#!/usr/bin/env python3

import argparse
import itertools
import logging

from multiprocessing import Pool
from collections import defaultdict
from operator import attrgetter
from typing import Collection, Iterator

from frippa.models import Protein, Cluster, Organism, Scaffold
from frippa.runners import hmmsearch, radar, signalp
from frippa.results import parse

from tqdm import tqdm

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s", datefmt="%H:%M:%S"
)

log = logging.getLogger(__name__)
log.setLevel(logging.INFO)


def find_DUF3328_proteins(organism, threads=1, evalue=0.1):
    """Find Proteins with DUF3328 hits in an Organism using hmmsearch."""

    log.debug("Launching hmmsearch")
    domtbl = hmmsearch(organism, threads=threads)

    log.debug("Parsing results")
    domtbl = parse(domtbl.split("\n"), "hmmsearch", evalue_cutoff=evalue)

    log.debug("Updating Protein objects")
    for scaffold in organism.records.values():
        for protein in scaffold.proteins:
            if protein.name in domtbl:
                protein.duf = domtbl[protein.name]


def translate_frame(
    sequence,
    frame: int,
    strand: int,
    transl_table: int,
    min_orf_length: int,
) -> Iterator[dict]:
    """Finds all ORFs for a given `frame` of `sequence`.

    Arguments:
        sequence: Nucleotide sequence
        frame: Translation frame
        strand: Sequence strand
        transl_table: Translation table
        min_orf_length: Minimum ORF length

    Yields:
        Dictionaries corresponding to predicted ORFs greater than
        `min_orf_length` residues long.
    """
    # Translate this frame ensuring multiples of 3
    frame_len = len(sequence) - frame
    frame_end = frame + frame_len - frame_len % 3
    trans = sequence[frame: frame_end].translate(transl_table)
    if not trans:
        return

    # Find start indices of ORFs via stop codons
    starts = [i + 1 for i, amino in enumerate(trans) if amino == "*"]
    starts = [0, *starts, len(trans)]

    # Iterate starts/ends
    for index in range(len(starts) - 1):
        aa_start, aa_end = starts[index: index + 2]
        orf_sequence = trans[aa_start: aa_end]

        # Want ORFs beginning with Methionine, ending *
        if not sequence[0] == "M":
            first_met = orf_sequence.find("M")
            if first_met == -1:
                continue
            orf_sequence = orf_sequence[first_met:]

        # Filter out any sequences that are too short
        if len(orf_sequence) < min_orf_length:
            continue

        # Correlate AA start/end to sequence, save protein
        yield dict(
            sequence=str(orf_sequence),
            start=frame + aa_start * 3,
            end=frame + aa_end * 3 + 3,
            strand=strand,
            frame=frame,
        )


def get_six_frame_translation(
    sequence,
    transl_table: int = 11,
    min_orf_length: int = 100
) -> Iterator[dict]:
    """Finds putative proteins from 6-frame translation of `cluster` sequence.

    Adapted from BioPython guide, section 20.1.13:
    http://biopython.org/DIST/docs/tutorial/Tutorial.html
    """
    seq_len = len(sequence)
    for strand in [+1, -1]:
        if strand == -1:
            sequence = sequence.reverse_complement()
        for frame in range(3):
            yield from translate_frame(
                sequence=sequence,
                frame=frame,
                strand=strand,
                transl_table=transl_table,
                min_orf_length=min_orf_length,
            )


def get_surrounding_proteins(middle, proteins, cutoff=10000):
    """Get Proteins surrounding an anchor Protein within some cutoff.

    Start with the anchor index, iterate either side until either out of cutoff range,
    or no more proteins left in that direction.
    """
    anchor = proteins[middle]
    neighbours = [anchor]

    # Stop flags and bounds check lambdas
    lflag, ltest = False, lambda x: x.end < anchor.start - cutoff
    uflag, utest = False, lambda x: x.start > anchor.end + cutoff

    count, total = 1, len(proteins)
    while not (lflag or uflag):
        lower, upper = middle - count, middle + count

        if not lflag and lower < 0:
            lflag = True
        if not uflag and upper >= total:
            uflag = True

        for index, flag, test in [(lower, lflag, ltest), (upper, uflag, utest)]:
            if flag:
                continue
            protein = proteins[index]
            if test(protein):
                flag = True
            if not flag:
                neighbours.append(protein)
        count += 1
    return sorted(neighbours, key=attrgetter("start"))


def initialise_clusters(organism, cutoff=10000):
    """Create Cluster objects for every DUF3328 protein."""
    for scaffold in organism.records.values():
        scaffold_len = len(scaffold.sequence)
        for index, protein in enumerate(scaffold.proteins):
            if not protein.duf:
                continue
            neighbours = get_surrounding_proteins(index, scaffold.proteins, cutoff=cutoff)
            start = neighbours[0].start
            end = neighbours[-1].end
            cluster = Cluster(
                proteins=neighbours,
                scaffold=scaffold.name,
                sequence=scaffold.sequence[start: end],
                start=start,
                end=end
            )
            organism.clusters.append(cluster)


def clusters_overlap(one, two):
    """Test if two clusters share the same precursor and should be merged."""
    one_peptides = set(one.proteins[i] for i in one.precursors)
    two_peptides = set(two.proteins[i] for i in two.precursors)
    return not one_peptides.isdisjoint(two_peptides)


def merge_overlapping_clusters(organism, overlap=0.6):
    """Merge Clusters in an Organism if they overlap more than some threshold."""
    index = 1
    while index < len(organism.clusters):
        one, two = organism.clusters[index - 1: index + 1]
        if clusters_overlap(one, two):
            log.debug("Merging %s and %s", one.location, two.location)
            one.add_proteins(two.proteins)
            del organism.clusters[index]
        else:
            index += 1
    for cluster in organism.clusters:
        cluster.proteins.sort(key=attrgetter("start"))


def find_signal_peptides(organism):
    proteins = list(
        itertools.chain.from_iterable(
            cluster.proteins for cluster in organism.clusters
        )
    )
    log.info("Running SignalP on %s proteins", len(proteins))
    output = parse(
        signalp(proteins).split("\n"),
        "signalp"
    )
    for protein in proteins:
        if protein.name not in output:
            continue
        protein.signalp = output[protein.name]


def find_repeats_in_sequence(
    protein,
    min_repeats=3,
    min_repeat_length=8,
    max_repeat_length=40,
    z_score_cutoff=10,
    cutsite_pct=0.6,
    repeat_similarity=0.8,
):
    repeats = []
    try:
        repeats = parse(
            radar(protein),
            "radar",
            min_repeats=min_repeats,
            min_repeat_length=min_repeat_length,
            max_repeat_length=max_repeat_length,
            z_score_cutoff=z_score_cutoff,
            cutsite_pct=cutsite_pct,
            repeat_similarity=repeat_similarity,
        )
    except UnicodeDecodeError:
        # bad continuation
        log.error("RADAR failed on %s", protein.name)
    return repeats


def find_repeats(
    organism: Organism,
    min_repeats: int=3,
    min_repeat_length: int=8,
    max_repeat_length: int=40,
    z_score_cutoff: int=10,
    cutsite_pct: float=0.6,
    repeat_similarity: float=0.8,
) -> None:
    """Detects repeat sequences in `organism` clusters using RADAR."""
    proteins = [
        protein
        for cluster in organism.clusters
        for protein in cluster.proteins
        if not protein.duf and protein.signalp
    ]
    log.info("Running RADAR on %s proteins", len(proteins))
    for protein in proteins:
        protein.repeats = find_repeats_in_sequence(
            protein,
            min_repeats=min_repeats,
            min_repeat_length=min_repeat_length,
            max_repeat_length=max_repeat_length,
            z_score_cutoff=z_score_cutoff,
            cutsite_pct=cutsite_pct,
            repeat_similarity=repeat_similarity,
        )


def find_unannotated_gaps(scaffold: Scaffold) -> Collection[tuple]:
    """Finds start and end coordinates of all inter-protein gaps in `scaffold`."""
    pairs = []
    total_proteins = len(scaffold.proteins)
    if total_proteins == 0:
        return pairs
    for index in range(total_proteins - 1):
        a, b = scaffold.proteins[index: index + 2]
        pair = (a.end + 1, b.start - 1)
        pairs.append(pair)
    scaffold_len = len(scaffold)
    if scaffold.proteins[0].start != 0:
        pairs.insert(0, (0, scaffold.proteins[0].start))
    if scaffold.proteins[-1].end != scaffold_len:
        pairs.append((scaffold.proteins[-1].end, scaffold_len))
    return pairs


def find_ORFs_in_scaffold(scaffold, transl_table: int=1, min_orf_length: int=100):
    """Find ORFs in unannotated gaps between proteins in `scaffold`."""
    count = 0
    for start, end in find_unannotated_gaps(scaffold):
        sequence = scaffold.sequence[start: end]
        for orf in get_six_frame_translation(
            sequence,
            transl_table=transl_table,
            min_orf_length=min_orf_length
        ):
            if not orf:
                continue
            orient = "F" if orf["strand"] == 1 else "R"
            protein = Protein(
                name=f"{scaffold.name}_ORF{count}_{orient}{orf['frame']}",
                scaffold=scaffold.name,
                start=start + orf["start"],
                end=start + orf["end"],
                strand=orf["strand"],
                sequence=orf["sequence"]
            )
            scaffold.proteins.append(protein)
            count += 1
    scaffold.proteins.sort(key=attrgetter("start"))
    return count


def find_ORFs_in_organism(organism: Organism, transl_table: int=1, min_orf_length: int=100) -> int:
    """Finds ORFs in scaffolds in `organism`, returning the count of added ORFs."""
    total_scaffolds = len(organism.records)
    pbar = tqdm(organism.records.values())
    count = 0
    for scaffold in pbar:
        pbar.set_description(scaffold.name, refresh=True)
        count += find_ORFs_in_scaffold(
            scaffold,
            transl_table=transl_table,
            min_orf_length=min_orf_length,
        )
    return count


def frippa(
    genbank,
    output=None,
    threads=1,
    evalue=0.1,
    neighbours=10000,
    overlap=0.8,
    min_repeats=3,
    min_repeat_length=8,
    max_repeat_length=40,
    cutsite_pct=0.8,
    report_all=False,
    z_score_cutoff=10,
    repeat_similarity=0.8,
    include_orfs=False,
    min_orf_length=100,
    transl_table=1
):
    """Main entry point to fRiPPa."""

    log.info("Parsing %s", genbank)
    organism = Organism.from_genbank(genbank)

    if include_orfs:
        log.info("Finding ORFs in genome >= %i AA", min_orf_length)
        count = find_ORFs_in_organism(
            organism,
            transl_table=transl_table,
            min_orf_length=min_orf_length,
        )
        log.info("Added %i ORFs", count)

    log.info("Searching for proteins with DUF3328 domains")
    find_DUF3328_proteins(organism, threads, evalue)

    log.info("Collecting proteins within %s bp of a DUF3328 hit", neighbours)
    initialise_clusters(organism, cutoff=neighbours)

    log.info("Searching for signal peptides in neighbour proteins")
    find_signal_peptides(organism)

    log.info("Searching for amino acid repeats in neighbour proteins")
    find_repeats(
        organism,
        min_repeats=min_repeats,
        min_repeat_length=min_repeat_length,
        max_repeat_length=max_repeat_length,
        z_score_cutoff=z_score_cutoff,
        cutsite_pct=cutsite_pct,
        repeat_similarity=repeat_similarity,
    )

    for cluster in organism.clusters:
        cluster.count()

    log.info("Finalising clusters")
    merge_overlapping_clusters(organism)

    if not report_all:
        organism.remove_invalid_clusters()

    if not organism.clusters:
        log.info("No clusters were predicted. Perhaps try different settings?")
    else:
        log.info("Generating summary")
        report = "\n\n".join(cluster.summary() for cluster in organism.clusters)
        if output:
            log.info("Writing to %s", output)
            with open(output, "w") as out:
                out.write(report)
        else:
            print(report)
    return organism


def get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument("genbank", help="Input GenBank file")
    parser.add_argument("-o", "--output", help="Path to write results to")
    parser.add_argument(
        "-t",
        "--threads",
        default=1,
        type=int,
        help="How many threads to use when running hmmsearch and SignalP (def. 1)",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        default=0.01,
        type=float,
        help="E-value cutoff to use when filtering hmmsearch results (def. 0.01).",
    )

    orfs = parser.add_argument_group("Open reading frame prediction")
    orfs.add_argument(
        "-orfs",
        "--include_orfs",
        action="store_true",
        help="Find ORFs between genes to include in analysis"
    )
    orfs.add_argument(
        "-mol",
        "--min_orf_length",
        type=int,
        default=100,
        help="Minimum ORF length (def. 100)."
    )
    orfs.add_argument(
        "-tt",
        "--transl_table",
        type=int,
        default=1,
        help="Translation table (def. 1)"
    )

    cluster = parser.add_argument_group("Cluster prediction settings")
    cluster.add_argument(
        "-n",
        "--neighbours",
        type=float,
        default=10000,
        help="Number of base pair either side of a DUF3328 hit to grab proteins (def. 10000).",
    )
    cluster.add_argument(
        "-v",
        "--overlap",
        type=float,
        default="0.6",
        help="Maximum percentage of overlapping proteins before clusters are merged (def. 0.6).",
    )
    cluster.add_argument(
        "-ra",
        "--report_all",
        action="store_true",
        help="Report all clusters with DUF3328 hits, even if no repeat predictions",
    )

    repeat = parser.add_argument_group("Repeat filtering")
    repeat.add_argument(
        "-cp",
        "--cutsite_pct",
        type=float,
        default=0.5,
        help="Threshold for total cut sites in individual repeats (def. 0.5).",
    )
    repeat.add_argument(
        "-mr",
        "--min_repeats",
        type=int,
        default=3,
        help="Minimum number of detected repeats (def. 3)."
    )
    repeat.add_argument(
        "-min_rl",
        "--min_repeat_length",
        type=int,
        default=8,
        help="Minimum length of a repeat sequence (def. 8)."
    )
    repeat.add_argument(
        "-max_rl",
        "--max_repeat_length",
        type=int,
        default=40,
        help="Maximum length of a repeat sequence (def. 40)."
        " Long sequences tend to be red herrings.",
    )
    repeat.add_argument(
        "-rs",
        "--repeat_similarity",
        type=float,
        default=0.5,
        help="Average (hamming) similarity threshold for repeat sequences (def. 0.5).",
    )
    repeat.add_argument(
        "-zc",
        "--z_score_cutoff",
        type=int,
        default=0,
        help="Minimum z-score for a repeat sequence (def. 0).",
    )
    return parser.parse_args()


def main():
    args = get_arguments()
    frippa(
        args.genbank,
        output=args.output,
        threads=args.threads,
        evalue=args.evalue,
        neighbours=args.neighbours,
        overlap=args.overlap,
        report_all=args.report_all,
        cutsite_pct=args.cutsite_pct,
        min_repeats=args.min_repeats,
        min_repeat_length=args.min_repeat_length,
        max_repeat_length=args.max_repeat_length,
        repeat_similarity=args.repeat_similarity,
        include_orfs=args.include_orfs,
        min_orf_length=args.min_orf_length,
        transl_table=args.transl_table
    )


if __name__ == "__main__":
    main()
