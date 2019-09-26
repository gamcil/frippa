#!/usr/bin/env python3

"""
Module for parsing results of each program.
"""

import re

from collections import defaultdict


class DomainHit:
    __slots__ = ("domain", "start", "end", "evalue")

    def __init__(self, domain, start, end, evalue):
        self.domain = domain
        self.start = int(start)
        self.end = int(end)
        self.evalue = float(evalue)

    def __repr__(self):
        return f"{self.domain} [{self.start}-{self.end}] {self.evalue}"


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
            domain=fields[2], start=fields[18], end=fields[19], evalue=fields[11]
        )

        if domain.evalue < evalue_cutoff:
            candidates[protein_id].append(domain)

    return candidates


class SignalPeptide:
    """The Repeat class represents individual repeat sequences from RADAR output.
    """

    __slots__ = ("start", "end", "site", "probability")

    def __init__(self, start, end, site, probability):
        self.start = int(start)
        self.end = int(end)
        self.site = site
        self.probability = float(probability)

    def __repr__(self):
        return "{}\t{}\t{}\t{}".format(
            self.start, self.end, self.site, self.probability
        )


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


class Match:
    """The Match class represents distinct blocks from RADAR output.
    """

    __slots__ = ("score", "length", "diagonal", "bw_from", "bw_to", "level", "repeats")

    def __init__(self, score, length, diagonal, bw_from, bw_to, level, repeats):
        self.score = float(score)
        self.length = int(length)
        self.diagonal = int(diagonal)
        self.bw_from = int(bw_from)
        self.bw_to = int(bw_to)
        self.level = int(level)
        self.repeats = repeats

    def __repr__(self):
        return "score: {}, length: {}, {} repeats".format(
            self.score, self.length, len(self.repeats)
        )


class Repeat:
    """The Repeat class represents individual repeat sequences from RADAR output.
    """

    __slots__ = ("start", "end", "score", "z_score", "sequence")

    def __init__(self, start, end, score, z_score, sequence):
        self.start = int(start)
        self.end = int(end)
        self.score = float(score)
        self.z_score = float(z_score)
        self.sequence = sequence

    def __repr__(self):
        return "{}\t{}\t{}\t{}\t{}".format(
            self.start, self.end, self.score, self.z_score, self.sequence
        )


def has_enough_cut_sites(repeats, threshold=0.5):
    """Check if number of cut sites in a list of repeats is over some threshold.
    """
    cut_sites = ["KR", "KK", "RR"]
    total = len(repeats)
    count = sum(site in repeat.sequence for site in cut_sites for repeat in repeats)
    return True if count / total >= threshold else False


def _radar(handle, max_repeat_length=40):
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

        if (  # Either no repeats, too long or not enough cut sites
            not repeats
            or any(len(repeat.sequence) > max_repeat_length for repeat in repeats)
            or not has_enough_cut_sites(repeats)
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


def parse(handle, program):
    """
    """
    parsers = {"hmmsearch": _hmmsearch, "radar": _radar, "signalp": _signalp}

    if program not in parsers:
        raise ValueError(f"Invalid program specified: {program}")

    return parsers[program](handle)
