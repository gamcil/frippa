#!/usr/bin/env python3

"""
"""

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

    def __init__(self, proteins):
        self.proteins = proteins

        self.duf = []
        self.precursors = []
        self.signalp = []

        self.count()

    def __repr__(self):
        return f"RiPP cluster [{len(self.proteins)} proteins, {self.start}-{self.end}]"

    @property
    def start(self):
        return self.proteins[0].start

    @property
    def end(self):
        return self.proteins[-1].end

    def count(self):
        """Iterate the Proteins in this Cluster and save indices of those with DUF3328
        hits, repeats (i.e. precursor peptides) and signal peptides (SignalP hits).
        """
        for index, protein in enumerate(self.proteins):
            if protein.duf:
                self.duf.append(index)
            if protein.repeats:
                self.precursor.append(index)
            if protein.signalp:
                self.signalp.append(index)

    def summarise(self, delimiter=","):
        """Build a Protein summary table for this Cluster.

        e.g.
            Index   Protein     DUF3328   Precursor  Signal Peptide
            1       GENE_0001   X         -          -
            2       GENE_0002   -         X          -
            3       GENE_0003   -         -          X
        """

        def check(index, array):
            return "X" if index in array else "-"

        def form_row(index, protein):
            checks = [
                check(index, array)
                for array in (self.duf, self.precursors, self.signalp)
            ]
            return delimiter.join([index, protein, *checks])

        header = ["Index", "Protein", "DUF3328", "Precursor", "Signal Peptide"]
        report = [form_row(index, protein) for index, protein in self.proteins]

        return "\n".join(header + report)

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
                Protein.from_SeqFeature(feature, record.seq)
                for feature in record.features
                if feature.type == "CDS"
            ]
            for record in SeqIO.parse(file, "genbank")
        }
        return cls(records)


class Protein:
    """Representation of a Protein.
    """

    def __init__(self, name, sequence, start, end, strand):
        self.name = name
        self.sequence = sequence
        self.start = start
        self.end = end
        self.strand = strand

        # Annotation status
        self.duf = False
        self.repeats = []
        self.signalp = []

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
            translation = feature.extract(record).translate()

        return cls(
            name=name[0],
            start=int(feature.location.start),
            end=int(feature.location.end),
            strand=int(feature.location.strand),
            sequence=translation,
        )
