#!/usr/bin/env python3

"""
"""

import json

from Bio import SeqIO


class Cluster:
    """Store proteins in a RiPP cluster, DUF3328 hits, precursor peptide.

    Parameters
    ----------
    _proteins : list
        Protein objects in this cluster.
    duf : list
        Indices in `proteins` of elements with DUF3328 hits.
    precursors : list
        Indices in `proteins` of elements with detected repeats (i.e. precursors)
    signalp : list
        Indicies in `proteins` of elements with signal peptides.
    """

    def __init__(self, proteins, scaffold):
        self.proteins = proteins
        self.scaffold = scaffold

        self.duf = set()
        self.precursors = set()
        self.signalp = set()

        self.count()

    def __str__(self):
        return "RiPP cluster [{} proteins, {}]".format(
            len(self.proteins), self.location
        )

    def to_json(self):
        pass

    @property
    def is_valid(self):
        return all([self.duf, self.precursors, self.signalp])

    @property
    def start(self):
        return self.proteins[0].start

    @property
    def end(self):
        return self.proteins[-1].end

    @property
    def location(self):
        return "{}:{}-{}".format(self.scaffold, self.start, self.end)

    def count(self):
        """Iterate the Proteins in this Cluster and save indices of those with DUF3328
        hits, repeats (i.e. precursor peptides) and signal peptides (SignalP hits).
        """
        for index, protein in enumerate(self.proteins):
            if protein.duf:
                self.duf.add(index)
            if protein.repeats:
                self.precursors.add(index)
            if protein.signalp:
                self.signalp.add(index)

    def summary(self, show_header=True):
        """Build a Protein summary table for this Cluster.

        e.g.
            Index   Protein     DUF3328   Precursor  Signal Peptide
            1       GENE_0001   X         -          -
            2       GENE_0002   -         X          -
            3       GENE_0003   -         -          X

        Pass 'hr' to `delimiter` for human readable format (i.e. spaced for pretty
        printing).
        """

        def check(index):
            return [
                "X" if index in array else "-"
                for array in (self.duf, self.precursors, self.signalp)
            ]

        def human_row(fields, protein_width):
            return "{: <8}{: <{}}{: <10}{: <12}{}".format(
                *fields[0:2], protein_width + 3, *fields[2:]
            )

        header = ["Index", "Protein", "DUF3328", "Precursor", "Signal Peptide"]
        report = [header] if show_header else []

        for index, protein in enumerate(self.proteins):
            row = [index, protein.name, *check(index)]
            report.append(row)

        protein_width = max(7, max(len(p.name) for p in self.proteins))
        report = [human_row(row, protein_width) for row in report]

        repeats = "\n\n".join(
            self.proteins[idx].name + "\n" + self.proteins[idx].summarise_repeats()
            for idx in self.precursors
        )

        if not repeats:
            repeats = "No repeats found..."

        separator = "-" * len(report[0])

        if show_header:
            report[0] = "{}\n{}\n{}".format(separator, report[0], separator)

        return "{}\n\n{}\n\n{}\n{}".format(self, repeats, "\n".join(report), separator)

    @property
    def fasta(self):
        """Return fasta record of all Proteins in the Cluster
        """
        return "\n".join([p.fasta for p in self.proteins])

    @property
    def proteins_with_sites(self):
        """Return proteins in this Cluster than contain Trypsin cut sites
        """
        return [
            protein
            for protein in self.proteins
            if "KK" in protein.translation
            or "KR" in protein.translation
            or "RR" in protein.translation
        ]

    def add_proteins(self, proteins):
        for protein in proteins:
            if protein in self.proteins:
                continue
            self.proteins.append(protein)
        self.count()


class Organism:
    """Representation of a genome file.
    """

    def __init__(self, records):
        self.records = records
        self.clusters = []

    def remove_invalid_clusters(self):
        self.clusters = [cluster for cluster in self.clusters if cluster.is_valid]

    @property
    def proteome(self):
        return "\n".join(
            protein.fasta for proteins in self.records.values() for protein in proteins
        )

    @classmethod
    def from_genbank(cls, file):
        """Parse a GenBank file using BioPython and create a new Organism.
        """
        records = {
            record.name: [
                Protein.from_SeqFeature(feature, record)
                for feature in record.features
                if feature.type == "CDS"
            ]
            for record in SeqIO.parse(file, "genbank")
        }
        return cls(records)


class Protein:
    """Representation of a Protein.
    """

    def __init__(self, name, sequence, scaffold, start, end, strand):
        self.name = name
        self.sequence = sequence
        self.scaffold = scaffold
        self.start = start
        self.end = end
        self.strand = strand

        # Annotation status
        self.duf = False
        self.repeats = []
        self.signalp = []

    def __str__(self):
        return "{} [{}:{}-{}({})]".format(
            self.name, self.scaffold, self.start, self.end, self.strand
        )

    def summarise_repeats(self):
        return "\n\n".join(
            f"{index: <5}{match.repeats[0].sequence}\n"
            + "\n".join(f"     {repeat.sequence}" for repeat in match.repeats[1:])
            for index, match in enumerate(self.repeats, 1)
        )

    @property
    def fasta(self):
        return f">{self.name}\n{self.sequence}"

    @classmethod
    def from_SeqFeature(cls, feature, record):
        """Create a Protein object from a BioPython SeqFeature object.

        SeqFeature passed to this function should be a 'CDS' type.
        """
        for identifier in ["protein_id", "locus_tag", "ID"]:
            name = feature.qualifiers.get(identifier, None)
            if name:
                break
        if not name:
            raise ValueError("Could not determine a sequence identifier")

        try:
            translation = feature.qualifiers["translation"][0]
        except KeyError:
            translation = feature.extract(record.seq).translate()

        return cls(
            name=name[0],
            scaffold=record.id,
            start=int(feature.location.start),
            end=int(feature.location.end),
            strand=int(feature.location.strand),
            sequence=translation,
        )


class DomainHit:
    __slots__ = ("domain", "start", "end", "evalue")

    def __init__(self, domain, start, end, evalue):
        self.domain = domain
        self.start = int(start)
        self.end = int(end)
        self.evalue = float(evalue)

    def __str__(self):
        return f"{self.domain} [{self.start}-{self.end}] {self.evalue}"


class SignalPeptide:
    """The Repeat class represents individual repeat sequences from RADAR output.
    """

    __slots__ = ("start", "end", "site", "probability")

    def __init__(self, start, end, site, probability):
        self.start = int(start)
        self.end = int(end)
        self.site = site
        self.probability = float(probability)

    def __str__(self):
        return "{}\t{}\t{}\t{}".format(
            self.start, self.end, self.site, self.probability
        )


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

    def __str__(self):
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

    def __str__(self):
        return "{}\t{}\t{}\t{}\t{}".format(
            self.start, self.end, self.score, self.z_score, self.sequence
        )
