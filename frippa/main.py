#!/usr/bin/env python3

import argparse
import itertools
import logging

from operator import attrgetter

from frippa.models import Cluster, Organism
from frippa.runners import hmmsearch, radar, signalp
from frippa.results import parse

logging.basicConfig(
    format="[%(asctime)s] %(levelname)s - %(message)s", datefmt="%H:%M:%S"
)

log = logging.getLogger(__name__)


def find_DUF3328_proteins(organism, threads=1, evalue=0.1):
    """Find Proteins with DUF3328 hits in an Organism using hmmsearch."""

    log.debug("Launching hmmsearch")
    domtbl = hmmsearch(organism, threads=threads)

    log.debug("Parsing results")
    domtbl = parse(domtbl.split("\n"), "hmmsearch", evalue_cutoff=evalue)

    log.debug("Updating Protein objects")
    for record in organism.records.values():
        for protein in record:
            if protein.name in domtbl:
                protein.duf = domtbl[protein.name]


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
    organism.clusters = [
        Cluster(
            get_surrounding_proteins(index, proteins, cutoff=cutoff),
            scaffold=scaffold
        )
        for scaffold, proteins in organism.records.items()
        for index, protein in enumerate(proteins)
        if protein.duf
    ]


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


def find_signal_peptides(organism):
    proteins = list(
        itertools.chain.from_iterable(
            cluster.proteins for cluster in organism.clusters
        )
    )
    log.debug("Running SignalP on %s proteins", len(proteins))
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
    organism,
    min_repeats=3,
    min_repeat_length=8,
    max_repeat_length=40,
    z_score_cutoff=10,
    cutsite_pct=0.6,
    repeat_similarity=0.8,
):
    proteins = [
        protein
        for cluster in organism.clusters
        for protein in cluster.proteins
        if not protein.duf and protein.signalp
    ]
    log.debug("Running RADAR on %s proteins", len(proteins))
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
):
    """Main entry point to fRiPPa."""

    log.info("Parsing %s", genbank)
    organism = Organism.from_genbank(genbank)

    log.info("Searching for proteins with DUF3328 domains")
    find_DUF3328_proteins(organism, threads, evalue)

    log.info("Collecting proteins within %s bp of a DUF3328 hit", neighbours)
    initialise_clusters(organism, cutoff=neighbours)

    log.info("Searching for signal peptides in neighbour proteins with SignalP")
    find_signal_peptides(organism)

    log.info("Searching for amino acid repeats in neighbour proteins with RADAR")
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
        help="How many threads to use when running hmmsearch and SignalP",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        default=0.01,
        type=float,
        help="E-value cutoff to use when filtering hmmsearch results",
    )

    cluster = parser.add_argument_group("Cluster prediction settings")
    cluster.add_argument(
        "-n",
        "--neighbours",
        type=float,
        default=10000,
        help="Number of base pair either side of a DUF3328 hit to grab proteins",
    )
    cluster.add_argument(
        "-v",
        "--overlap",
        type=float,
        default="0.6",
        help="Maximum percentage of overlapping proteins before clusters are merged",
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
        help="Threshold for total cut sites in individual repeats",
    )
    repeat.add_argument(
        "-mr",
        "--min_repeats",
        type=int,
        default=3,
        help="Minimum number of detected repeats"
    )
    repeat.add_argument(
        "-min_rl",
        "--min_repeat_length",
        type=int,
        default=8,
        help="Minimum length of a repeat sequence."
    )
    repeat.add_argument(
        "-max_rl",
        "--max_repeat_length",
        type=int,
        default=40,
        help="Maximum length of a repeat sequence."
        " Long sequences tend to be red herrings",
    )
    repeat.add_argument(
        "-rs",
        "--repeat_similarity",
        type=float,
        default=0.8,
        help="Average (hamming) similarity threshold for repeat sequences",
    )
    repeat.add_argument(
        "-zc",
        "--z_score_cutoff",
        type=int,
        default=0,
        help="Minimum z-score for a repeat sequence",
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
    )


if __name__ == "__main__":
    main()
