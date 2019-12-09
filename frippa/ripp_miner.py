#!/usr/bin/env python3

"""
Identify potential RiPP clusters in funcoDB genomes.

1. Search for proteins with DUF3328 Pfam domain,
   either directly by querying pre-annotated genomes,
   or making single profile HMM database and running
   hmmer on proteomes

2. Take surrounding +- 10 proteins (user threshold)

3. Filter based on presence of signal peptide and
   trypsin cut site, i.e. SignalP and KR/KK

4. Search for repeats in remaining protein sequences
   using RADAR

5. Return TSV:
       a) Protein with DUF3328
       b) Nearby proteins that match criteria, distance
          from DUF3328 protein, functional annotations
       c) Annotations of other genes in proximity
"""

import sys
import logging
import sqlite3
import subprocess

from pathlib import Path
from multiprocessing import Pool
from tempfile import NamedTemporaryFile

from funcodb import config, utils
from funcodb.classes import Cluster


class RiPPCluster(Cluster):
    """Store proteins in a RiPP cluster, DUF3328 hits, precursor peptide
    """
    def __init__(self, duf_hit):
        # Even though it takes no args to instantiate a Cluster, still
        # want to have basic fields defined
        super().__init__()

        self.duf_hits = [duf_hit]   # IDs of Proteins with DUF3328 hits
        self.precursors = {}        # Proteins with detected repeats
        self.signalp_hits = []      # Proteins with signal peptides

    @property
    def valid_cluster(self):
        """Simple flag that returns true if this Cluster has all necessary
        components of the RiPP cluster (each list is non-empty)
        """
        return self.duf_hits and self.precursors and self.signalp_hits

    @property
    def fasta(self):
        """Return fasta record of all Proteins in the Cluster
        """
        return '\n'.join([p.fasta for p in self.proteins])

    @property
    def proteins_with_sites(self):
        """Return proteins in this Cluster than contain Trypsin cut sites
        """
        return [protein for protein in self.proteins if
                'KK' in protein.translation or
                'KR' in protein.translation or
                'RR' in protein.translation]

    # TODO: just make base assimilate() method take a list of attribute
    #       lists to merge isntead of overriding
    def assimilate(self, cluster_two):
        """Merge two RiPP clusters.

        Override of Cluster.assimilate(), only difference being we need to
        combine duf_hits, not core_genes
        """
        assert isinstance(cluster_two, Cluster)
        assert cluster_two.proteins is not None

        # Swallow new proteins and set new cluster boundaries
        for protein in cluster_two.proteins:
            if protein not in self.proteins:
                # Update the Cluster for this Protein
                protein.clusters.remove(cluster_two)
                protein.clusters.append(self)
                # Add to the Protein list
                self.proteins.append(protein)
        self.compute_boundaries_from_proteins()

        # Save DUF hits
        for hit in cluster_two.duf_hits:
            if hit not in self.duf_hits:
                self.duf_hits.append(hit)


def run_hmmsearch(hmmdb, organism, threshold=1.0):
    """Search for DUF3328 profile hits in proteome and return hits

    Want to run this in parallel using Pool, so limit --cpu to 1.
    """
    with NamedTemporaryFile() as proteome, NamedTemporaryFile() as table:
        # Write proteome of current organism to temp file
        proteome.write(organism.proteome_fasta.encode())
        command = [config.parser['Programs']['hmmsearch'],
                   '--cpu', '1',
                   '--domtblout', table.name,
                   hmmdb, proteome.name]

        # Run hmmsearch
        subprocess.run(command, stdout=subprocess.DEVNULL)

        # Parse hits table
        candidates = []
        for line in table:
            # Decode bytes, skip headers
            try:
                utf_line = line.decode()
            except UnicodeDecodeError:
                print('unicode error', line)
                continue
            if utf_line.startswith('#'):
                continue

            # Split line, grab protein ID and evalue
            # .split() will split on any whitespace and disregard ''
            fields = utf_line.strip().split()
            protein_id = fields[0]
            evalue = float(fields[12])

            # If evalue is lower than the supplied threshold, save
            if evalue < threshold:
                if protein_id not in candidates:
                    candidates.append(protein_id)

    return candidates


def predict_clusters(args):
    """Initial prediction of RiPP clusters in an Organism object.

        1. Run hmmsearch on the proteome, save DUF3328 hits
        2. Make Cluster objects, populate with neighbouring proteins
        3. Combine overlapping BGCs, save final list to Organism

    In multiprocessing.Pool, instances are pickled before being sent to
    the function, and thus the function just edits a copy of the object.

    Ergo, have to return the Organism object, can't just mutate it.
    """
    # Unpack the input tuple, run hmmsearch first
    hmmdb, organism, e_cutoff, neighbours, max_shared = args
    duf_hits = run_hmmsearch(hmmdb,     # Path to HMMer database
                             organism,  # Organism object
                             e_cutoff)  # DUF hit e-value cutoff
    for hit in duf_hits:
        # Grab Protein, Scaffold, and index of Proteins on the Scaffold
        protein = organism.protein_dict[hit]
        scaffold = protein.scaffold
        index, proteins = scaffold.ordered_proteins
        hit_index = index[hit]

        # Create Cluster object for each DUF hit, add neighbour Proteins
        cluster = RiPPCluster(hit)
        up_border = False
        down_border = False
        counter = 0
        while counter < neighbours:
            # If we've hit scaffold edge (no more proteins), break
            if up_border and down_border:
                break

            # Upstream
            if not up_border:
                try:
                    protein = proteins[hit_index + counter]
                    protein.clusters.append(cluster)
                    cluster.proteins.append(protein)
                except IndexError:
                    up_border = True

            # Downstream
            if not down_border:
                try:
                    protein = proteins[hit_index - counter]
                    protein.clusters.append(cluster)
                    cluster.proteins.append(protein)
                except IndexError:
                    down_border = True

            counter += 1

        # Add this Cluster to the Organism
        cluster.compute_boundaries_from_proteins()
        scaffold.clusters.append(cluster)

    # Scan scaffolds for overlapping Clusters and combine them
    for scaffold in organism.scaffolds:
        scaffold.merge_overlapping_clusters(max_shared)

        # De-reference irrelevant Protein instances
        scaffold.proteins = scaffold.cluster_proteins

    return organism


def run_signalp(log, proteins) -> list:
    """Run SignalP on a list of Proteins, return list of matches.

    Pool all proteins, analyse in one run, then parse results and map
    back to their DUF3328 hits.

    SignalP 5.0 is multithreaded but does not take a threads/cpus arg,
    only -batch, which controls how many sequences it will try and analyse
    simultaneously. Pre-SignalP protein pool for 6 genomes ~1600, so
    will likely have to think about changing this arg when scaling.
    """
    # Collate protein sequences into one FASTA for writing
    sequences = ''
    # SignalP changes | to _ in output so have to catch..
    flag = {}
    for protein in proteins:
        if '|' in protein.protein_id:
            changed = protein.protein_id.replace('|', '_')
            flag[changed] = protein.protein_id

        sequences += protein.fasta

    with NamedTemporaryFile() as fasta:
        # Write proteins to temp file and run SignalP
        fasta.write(sequences.encode())
        command = [config.parser['Programs']['signalp'],
                   '-fasta', fasta.name,
                   '-format', 'short',
                   '-prefix', 'funco']
        log.debug(f'Running SignalP 5.0 with {" ".join(command)}')
        print(f'SignalP on {len(proteins)} proteins')
        subprocess.run(command, stdout=subprocess.DEVNULL)

    # Parse the results, 
    output = Path(f'funco_summary.signalp5')
    candidates = []
    with open(output, 'r') as results:
        for line in results:
            if line.startswith('#') or 'OTHER' in line:
                continue
            # Change back altered protein IDs, add to list
            protein_id = line.strip().split('\t', 1)[0]
            if protein_id in flag:
                protein_id = flag[protein_id]

            candidates.append(protein_id)

    # Delete results file
    output.unlink()
    return(candidates)


def run_lfasta(sequence, matrix):
    """Run lfasta on protein sequence with given scoring matrix.

    RADAR combines output of BLOSUM50 and PAM250 matrices for
    repeat detection.

    This function must be called under a NamedTemporaryFile()
    context manager.

    Args:
        sequence (NamedTemporaryFile): path of tempfile containing
                                       amino acid sequence
        matrix (str): which scoring matrix to use
    Returns:
        (str): stdout from subprocess call
    """
    assert matrix in ['BLOSUM50', 'PAM250']

    # lfasta sequence sequence
    command = [config.parser['Programs']['lfasta'],
               sequence, sequence]

    # Append extra commands to use PAM250 matrix
    if matrix == 'PAM250':
        command.append('-f -12 -g -2 -s 250')

    # Run lfasta, capture and return stdout
    return subprocess.check_output(command)


def run_radar(input_tuple):
    """Run RADAR on list of proteins in single sequence mode

    -P  filename of multiple alignment [provide -S file here too]
    -Q  filename of lfasta output
    -R  filename of sequence
    -S  filename of sequence lfasta file [optional]
    -V  level of verbosity

    Args:
        proteins (list): Protein instances being analysed
    Returns:
        candidates (tuple): Proteins with repeats
    """
    # Unpack map tuple
    protein, max_length, min_cut_pct = input_tuple

    # Set up some tmpfiles for lfasta/RADAR to read
    with NamedTemporaryFile() as blosum, \
            NamedTemporaryFile() as pam, \
            NamedTemporaryFile() as sequence:
        # Protein sequence
        sequence.write(protein.translation.encode())
        sequence.flush()

        # Run lfasta with both BLOSUM50 and PAM250 matrices
        blosum.write(run_lfasta(sequence.name, 'BLOSUM50'))
        pam.write(run_lfasta(sequence.name, 'PAM250'))

        # Run RADAR and capture stdout
        command = [config.parser['Programs']['radar'],
                   '-P', sequence.name,
                   '-R', sequence.name,
                   '-Q', blosum.name,
                   '-S', pam.name]

        # TODO: error check subprocess calls; lfasta too
        try:
            results = subprocess.check_output(command).decode()
        except UnicodeDecodeError:
            print('unicode error in radar')
            return None

    # Split raw results on 75 dashes (RADAR output); easiest
    # way to isolate repeat blocks
    bad_outcomes = ['No repeats found',
                    'Not enough dots given']
    repeats = []
    for block in results.split(75*'-'):
        # Failure/skip cases
        if any(error in block for error in bad_outcomes):
            break
        elif not block or block == '\n' or \
                'No. of Repeats' in block:
            continue
        # Split on newline, then tab, to separate repeats
        repeats.append([rep.split('\t')[-1] for rep in
                       block.split('\n')[1:-1]])

    if not repeats:
        return None

    # Final check for length and cut sites
    cut_sites = ['KK', 'KR', 'RR']
    cuts = 0
    filtered = []
    for group in repeats:
        for repeat in group:
            # Skip long predicted repeats; likely wrong
            # TODO: this isn't working for some reason
            if len(repeat) > max_length:
                break

            # Tally up cuts
            # TODO: make this check >= 1 p/repeat, rather
            #       than just general percentage
            #       ie 3 reps -> [1, 0, 1] -> 2/3
            if any(site in repeat for site in cut_sites):
                cuts += 1

        # Save if cut sites appear in at least <threshold>
        if len(group) == 0:
            continue
        if cuts / len(group) >= min_cut_pct:
            filtered.append(group)

    return protein, filtered


def format_repeat(protein_id, repeats):
    """Format repeats nicely for printing.

        Protein_ID    Repeat1
                      Repeat2
                      Repeat3
    """
    margin = len(protein_id) + 4
    output = ''
    for i, repeat in enumerate(repeats):
        if i == 0:
            output = f'{protein_id}    {repeat}\n'
        else:
            space = ' ' * margin
            output += f'{space}{repeat}\n'
    return output


def main():
    """
    Run RiPP-Miner.

    Args:
        threads (int): number of CPUs to use
        e_cutoff (int): e-value threshold for hmmsearch
        n_cutoff (int): number of proteins to analyse either side of DUF hit
        max_length (int): maximum length of a detected repeat to keep
        max_shared (int): maximum number of proteins BGCs can share
                          without merging
        min_cut_pct (float): minimum number of individual repeats that must
                             contain a trypsin cut site
    """
    # TODO: replace with argparse
    threads = 1
    e_cutoff = 1.0
    n_cutoff = 15
    max_length = 40
    max_shared = 10
    min_cut_pct = 0.5

    # Set up logger and funcoDB path
    log = logging.getLogger(__name__)
    funcodb_path = config.parser['Database']['funcodb']

    # Path to the hmmer database containing the DUF3328 domain
    base_dir = utils.get_project_root()
    hmmdb = base_dir / 'funcodb' / 'ripp_miner' / 'DUF3328.hmm'

    # Get list of all organisms in funcoDB, then extract proteomes
    with sqlite3.connect(funcodb_path) as funco:
        print('Connected to funcoDB.')
        cursor = funco.cursor()

        # Get list of organism row IDs, then populate Organism instances
        print('Extracting organisms...')
        organism_ids = utils.list_organism_ids(cursor)
        organisms = {}  # Dict to store all Organism objects
        tuples = []     # Argument tuples for pool.imap()
        for oid in organism_ids[40:60]:
            # TODO: collapse this into Pool, as SQLite allows for many
            #       concurrent READERS, just not WRITERS
            #       then tasks = len(organism_ids)
            organisms[oid] = utils.organism_from_funco(cursor, oid)
            tuples.append((hmmdb, organisms[oid],
                           e_cutoff, n_cutoff,
                           max_shared))

        # Set up multiprocessing.Pool() to analyse Organisms in parallel
        with Pool(threads) as pool:
            tasks = len(tuples)
            # predict_clusters() MUST return Organism object, as objects
            #   sent into multiprocessing.Pool are pickled such that the
            #   daughter process essentially just edits a copy of the parent
            # imap is used instead of imap_unordered, as we want to preserve
            #   the order of organisms so we can re-assign them
            print(f'Predicting clusters in {oid}')
            for i, organism in enumerate(pool.imap(predict_clusters,
                                                   tuples), 1):
                organisms[i] = organism
                sys.stdout.write('Scanning organisms for clusters: '
                                 f'{i / tasks:.2%}\r')
                sys.stdout.flush()
            sys.stdout.write('\n\n')

        # Iterate organisms, and:
        # 1. Tally DUF3328 hits and clusters formed per genome, report
        # 2. Pool together all neighbour proteins for single SignalP run
        # 3. Populate dictionary of all protein instances for lookups
        print('ID   DUF3328   BGCs\tOrganism')
        print('-' * 55)
        neighbours = set()  # Avoid duplicates from Proteins in >1 BGC
        proteins = {}
        for oid in organisms:
            # Should now be a reduced dict, since we pruned scaffolds
            # just prior in predict_clusters()
            proteins.update(organisms[oid].protein_dict)
            clusters = 0
            duf_hits = 0
            for scaffold in organisms[oid].scaffolds:
                if scaffold.clusters:
                    clusters += len(scaffold.clusters)
                    for cluster in scaffold.clusters:
                        duf_hits += len(cluster.duf_hits)
                        neighbours.update(cluster.proteins)
            print(f'{oid}\t{duf_hits}\t{clusters}\t{organisms[oid].full_name}')

        # Analyse filtered proteins with SignalP
        # TODO: for in-house genomes, should already have SignalP
        #       annotations; for NCBI no and for JGI not yet
        log.info('Running SignalP on these proteins...')
        signalp = run_signalp(log, neighbours)
        print(f'{len(signalp)} proteins left from {len(neighbours)}.')

        # Run RADAR on any remaining proteins
        if signalp:
            log.info(f'SignalP found {len(signalp)} proteins with '
                     'signal peptides. Now running RADAR to detect '
                     'amino acid repeats...')

            # Form tuples for imap() and add proteins with signal peptides
            # to the Cluster list
            tuples = []
            for pid in signalp:
                protein = proteins[pid]
                tuples.append((protein, max_length, min_cut_pct))
                # TODO: fix this
                # protein.cluster.signalp_hits.append(protein)

            # tuples = [(proteins[pid], max_length, min_cut_pct) for
            #          pid in signalp]
            tasks = len(tuples)
            repeats = {}

            # Run RADAR, report progress
            with Pool(threads) as pool:
                for i, result in enumerate(pool.imap(run_radar, tuples), 1):
                    # Enumerate packages multiple returned values into tuple
                    # Unpack, then add to a dictionary for mapping outside Pool
                    if result:
                        protein, groups = result
                        repeats[protein.protein_id] = groups

                    sys.stdout.write('Scanning proteins for repeats: '
                                     f'{i / tasks:.2%}\r')
                    sys.stdout.flush()
                sys.stdout.write('\n\n')
        else:
            log.warning('SignalP returned no results!')
            raise SystemExit

        # Add repeats to corresponding clusters. Have to do outside Pool
        for pid in repeats:
            protein = proteins[pid]
            for cluster in protein.clusters:
                cluster.precursors[pid] = repeats[pid]

        # Print out
        with open('ripps.txt', 'w') as ripp:
            for oid in organisms:
                print(f'--- {oid}. {organisms[oid].full_name} ---')
                ripp.write(f'--- {oid}. {organisms[oid].full_name} ---\n')
                cluster_number = 1
                for scaffold in organisms[oid].scaffolds:
                    if not scaffold.clusters:
                        continue

                    print(f'\nScaffold: {scaffold.accession}')
                    ripp.write(f'\nScaffold: {scaffold.accession}\n')
                    for cluster in scaffold.clusters:
                        duf_hits = ', '.join(cluster.duf_hits)

                        print(f'Cluster {cluster_number}, '
                              f'DUF3328 hits: {duf_hits}')
                        ripp.write(f'Cluster {cluster_number}, '
                                   f'DUF3328 hits: {duf_hits}\n')
                        for pid in cluster.precursors:
                            for group in cluster.precursors[pid]:
                                print(format_repeat(pid, group))
                                ripp.write(format_repeat(pid, group) + '\n')

                        cluster_number += 1

def get_args():
    return

if __name__ == '__main__':
    main()
