#!/usr/bin/env python3


"""
This module provides the command line interface, as well as main routines for fRiPPa.
"""

import argparse
import itertools
import logging

from operator import attrgetter

from models import Cluster, Organism
from runners import hmmsearch, radar, signalp
from results import parse

logging.basicConfig(level=logging.DEBUG)

log = logging.getLogger(__name__)


def find_DUF3328_proteins(organism):
    """Find Proteins with DUF3328 hits in an Organism using hmmsearch."""

    log.debug("Launching hmmsearch")
    domtbl = hmmsearch(organism)

    log.debug("Parsing results")
    domtbl = parse(domtbl.split("\n"), "hmmsearch")

    log.debug("Updating Protein objects")
    for record in organism.records.values():
        for protein in record:
            if protein.name in domtbl:
                protein.duf = domtbl[protein.name]


def get_surrounding_proteins(middle, proteins, cutoff=10000):
    """Get Proteins surrounding an anchor Protein within some cutoff.
    """
    anchor = proteins[middle]
    neighbours = [anchor]

    lflag, ltest = False, lambda x: x.end < anchor.start - cutoff
    uflag, utest = False, lambda x: x.start > anchor.end + cutoff

    count, total = 1, len(proteins)
    while not (lflag or uflag):
        lower, upper = middle - count, middle + count

        if lower < 0:
            lflag = True
        if upper >= total:
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


def calculate_overlap(one, two):
    """Calculate overlap between two lists with respect to the first list."""
    overlap = set(one).intersection(two)
    return len(overlap) / len(one)


def initialise_clusters(organism, neighbour=10000):
    """Create Cluster objects for every DUF3328 protein."""
    for record in organism.records.values():
        for index, protein in enumerate(record):
            if not protein.duf:
                continue
            nearby = get_surrounding_proteins(index, record, cutoff=neighbour)
            cluster = Cluster(nearby)
            organism.clusters.append(cluster)


def filter_overlapping_clusters(organism, overlap=0.6):
    """Merge Clusters in an Organism if they overlap more than some threshold."""
    index = 1
    while index < len(organism.clusters):
        one, two = organism.clusters[index - 1 : index + 1]
        if calculate_overlap(one.proteins, two.proteins) > overlap:
            one.add_proteins(two.proteins)
            del organism.clusters[index]
        else:
            index += 1


def find_signal_peptides(organism):
    proteins = list(
        itertools.chain.from_iterable(cluster.proteins for cluster in organism.clusters)
    )

    log.debug("Running SignalP on %s proteins", len(proteins))
    output = parse(signalp(proteins).split("\n"), "signalp")

    for protein in proteins:
        if protein.name not in output:
            continue
        protein.signalp = output[protein.name]


def find_repeats(organism):
    for cluster in organism.clusters:
        for protein in cluster.proteins:
            if protein.duf:
                continue
            protein.repeats = parse(radar(protein), "radar")


def frippa(
    genbank,
    threads=1,
    evalue=0.1,
    neighbours=10,
    overlap=0.8,
    repeat_max=40,
    cutsite_pct=0.8,
):

    log.info("Parsing GenBank file")
    organism = Organism.from_genbank(genbank)

    log.info("Searching for proteins with DUF3328 domains")
    find_DUF3328_proteins(organism)

    log.info("Collecting proteins within %s proteins of a DUF3328 hit", neighbours)
    initialise_clusters(organism)
    filter_overlapping_clusters(organism)

    log.info("Searching for signal peptides in neighbour proteins with SignalP")
    find_signal_peptides(organism)

    log.info("Searching for amino acid repeats in neighbour proteins with RADAR")
    find_repeats(organism)

    log.info("Finalising potential RiPP clusters")

    return organism


def get_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help="How many threads to use when running hmmsearch and SignalP",
    )
    parser.add_argument(
        "-e",
        "--evalue",
        type=float,
        help="E-value cutoff to use when filtering hmmsearch results",
    )

    cluster = parser.add_argument_group("Cluster prediction settings")
    cluster.add_argument(
        "-n",
        "--neighbours",
        type=float,
        help="Number of proteins to grab either side of a DUF3328 hit",
    )
    cluster.add_argument(
        "-o",
        "--overlap",
        type=float,
        help="Maximum percentage of overlapping proteins before clusters are merged",
    )

    repeat = parser.add_argument_group("Repeat filtering")
    repeat.add_argument(
        "-r", "--repeat_max", type=int, help="Maximum length of a detected repeat"
    )
    repeat.add_argument(
        "-c",
        "--cutsite_pct",
        type=float,
        help="Threshold for total cut sites in individual repeats",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = get_arguments()
    frippa(**args)
