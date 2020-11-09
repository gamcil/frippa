#!/usr/bin/env python3


from Bio import SeqIO

from frippa import formatters


class Cluster:
    """Stores proteins in a RiPP cluster, DUF3328 hits, precursor peptide.

    Parameters:
        _proteins (list): Protein objects in this cluster.
        duf (list): Indices in proteins of elements with DUF3328 hits.
        precursors (list): Indices in proteins of elements with detected repeats (i.e. precursors)
        signalp (list): Indices in proteins of elements with signal peptides.
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
        repeats_fine = [repeat in self.signalp for repeat in self.precursors]
        return all([self.duf, repeats_fine])

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
        return formatters.format_cluster(self, show_headers=show_header)

    @property
    def fasta(self):
        """Builds FASTA record of all Proteins in the Cluster."""
        return "\n".join([p.fasta for p in self.proteins])

    @property
    def proteins_with_sites(self):
        """Returns proteins in this Cluster than contain Trypsin cut sites."""
        sites = ["KK", "KR", "RR"]
        return [
            protein
            for protein in self.proteins
            if any(site in protein.translation for site in sites)
        ]

    def add_proteins(self, proteins):
        for protein in proteins:
            if protein in self.proteins:
                continue
            self.proteins.append(protein)
        self.count()


class Organism:
    """An organism."""

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
    """A Protein sequence."""

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
    """A HMM profile domain hit."""

    __slots__ = ("domain", "start", "end", "evalue")

    def __init__(self, domain, start, end, evalue):
        self.domain = domain
        self.start = int(start)
        self.end = int(end)
        self.evalue = float(evalue)

    def __str__(self):
        return f"{self.domain} [{self.start}-{self.end}] {self.evalue}"


class SignalPeptide:
    """Signal peptide hits from SignalP."""

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
    """Distinct blocks from RADAR output."""

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
    """Individual repeat sequences from RADAR output."""

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
