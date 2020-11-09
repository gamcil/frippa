#!/usr/bin/env python3

"""
Module for parsing results of each program.
"""

import re

from collections import defaultdict

from frippa.models import DomainHit, SignalPeptide, Match, Repeat


def _hmmsearch(handle, evalue_cutoff=0.1):
    """Parse hmmsearch --domtblout file.

    This function will save any protein satisfying the evalue cutoff, as well as the
    name of the domain that was hit.
    """
    candidates = defaultdict(list)

    for line in handle:
        try:
            line = line.decode()
        except AttributeError:
            pass
        if not line or line.startswith("#") or line.isspace():
            continue

        protein_id, *fields = line.strip().split()

        domain = DomainHit(
            domain=fields[2],
            start=int(fields[18]),
            end=int(fields[19]),
            evalue=float(fields[11])
        )

        if domain.evalue < evalue_cutoff:
            candidates[protein_id].append(domain)

    return candidates


def _signalp(handle):
    """Parse output from SignalP 5.0.
    """
    hits = defaultdict(list)

    pattern = re.compile(
        r"CS pos: (?P<start>\d+)-(?P<end>\d+)\. (?P<site>\w+-\w+)\. Pr: (?P<prob>\d\.\d+)"
    )

    for line in handle:
        if not line or line.startswith("#") or "OTHER" in line:
            continue

        protein_id, _, _, _, information = line.strip().split("\t")

        search = pattern.search(information)

        peptide = SignalPeptide(
            search.group("start"),
            search.group("end"),
            search.group("site"),
            search.group("prob"),
        )

        hits[protein_id].append(peptide)

    return hits


def has_enough_cut_sites(repeats, threshold=0.5):
    """Check if number of cut sites in a list of repeats is over some threshold.
    """
    cut_sites = ["KR", "KK", "RR"]
    total = len(repeats)
    count = sum(site in repeat.sequence for site in cut_sites for repeat in repeats)
    return True if count / total >= threshold else False


def repeats_are_similar(repeats, threshold=0.5):
    """Compute similary between a collection of repeats.

    Computes hamming distance between each pair of repeat sequences, then divides by the
    number of repeats supplied (-1). Returns True if this score is over the specified
    threshold.
    """
    score, total, length = 0, len(repeats), len(repeats[0].sequence)
    for index in range(1, total):
        one, two = repeats[index - 1 : index + 1]
        score += sum(a == b for a, b in zip(one.sequence, two.sequence)) / length
    return score / (total - 1) >= threshold


def _radar(
    handle,
    min_repeats=3,
    min_repeat_length=8,
    max_repeat_length=40,
    z_score_cutoff=10,
    cutsite_pct=0.5,
    repeat_similarity=0.8,
):
    """Parse RADAR results.

    Uses a multiline regex to match result blocks, i.e.

    <ex>

    And then another regex to parse the actual repeat block, i.e.

    <ex>
    """

    stats_pattern = re.compile(
        r"-{75}\n"
        r"No\. of Repeats.+?\n"
        r".+?(?P<reps>\d+)\|"
        r".+?(?P<score>\d+.\d+)\|"
        r".+?(?P<length>\d+)\|"
        r".+?(?P<diag>\d+)\|"
        r".+?(?P<bw_from>\d+)\|"
        r".+?(?P<bw_to>\d+)\|"
        r".+?(?P<level>\d+).*?\n"
        r"-{75}\n"
        r"(?P<block>.+?)"
        r"\n-{75}\n",
        re.DOTALL,
    )

    block_pattern = re.compile(
        r"\s*?(?P<start>\d+)-"
        r"\s+(?P<end>\d+) "
        r"\((?P<score>\d+.\d+)/"
        r"(?P<zscore>\d+.\d+)\)\t"
        r"(?P<sequence>.+?)(?:\n|\Z)"  # non-capturing alternation group to use \z
    )

    matches = []
    for match in stats_pattern.finditer(handle):
        repeats = [
            Repeat(
                repeat.group("start"),
                repeat.group("end"),
                repeat.group("score"),
                repeat.group("zscore"),
                repeat.group("sequence"),
            )
            for repeat in block_pattern.finditer(match.group("block"))
        ]

        if (  # Either no repeats, too long, score too low, or not enough cut sites
            len(repeats) < min_repeats
            or any(
                len(repeat.sequence) < min_repeat_length
                or len(repeat.sequence) > max_repeat_length
                or repeat.z_score < z_score_cutoff
                for repeat in repeats
            )
            or not has_enough_cut_sites(repeats, threshold=cutsite_pct)
            or not repeats_are_similar(repeats, threshold=repeat_similarity)
        ):
            continue

        match = Match(
            match.group("score"),
            match.group("length"),
            match.group("diag"),
            match.group("bw_from"),
            match.group("bw_to"),
            match.group("level"),
            repeats,
        )
        matches.append(match)
    return matches


def parse(handle, program, **kwargs):
    parsers = {"hmmsearch": _hmmsearch, "radar": _radar, "signalp": _signalp}

    if program not in parsers:
        raise ValueError(f"Invalid program specified: {program}")

    return parsers[program](handle, **kwargs)
