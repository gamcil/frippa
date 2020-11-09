#!/usr/bin/env python3

"""
Module for handling subprocess calls for hmmsearch, lfasta, RADAR and SignalP.
"""

import os
import shutil
import subprocess

from pathlib import Path
from tempfile import NamedTemporaryFile as NTF


ROOT_DIR = Path(__name__).resolve().parent
DATA_DIR = ROOT_DIR / "data"


def get_path(*aliases):
    """Get path to a Program."""

    for alias in aliases:
        path = shutil.which(alias)
        if path:
            return path
    raise FileNotFoundError(f"Could not find {aliases} on system $PATH")


def make_fasta(proteins):
    return "\n".join(protein.fasta for protein in proteins)


def make_command(params):
    """Flatten parameter dict to list for command."""
    command = [params.pop("path")]
    last = params.pop("last", None)
    for key, value in params.items():
        if isinstance(value, list):
            command.extend([key, *value])
        else:
            command.extend([key, value])
    if isinstance(last, list):
        command.extend(last)
    elif isinstance(last, str):
        command.append(last)
    return command


def hmmsearch(organism, threads=1):
    """Search for DUF3328 profile hits in proteome and return hits."""

    params = {
        "path": get_path("hmmsearch"),
        "--cpu": str(threads),
        "-o": "/dev/null",
        "--domtblout": "/dev/stdout",
        "last": [DATA_DIR / "DUF3328.hmm", "-"],
    }

    return subprocess.run(
        make_command(params),
        input=organism.proteome,
        stdout=subprocess.PIPE,
        universal_newlines=True,
    ).stdout


def signalp(proteins, batch=1000):
    """Run SignalP on a list of Proteins, return list of matches.

    Pool all proteins, analyse in one run, then parse results and map
    back to their DUF3328 hits.

    SignalP 5.0 is multithreaded but does not take a threads/cpus arg,
    only -batch, which controls how many sequences it will try and analyse
    simultaneously. Pre-SignalP protein pool for 6 genomes ~1600, so
    will likely have to think about changing this arg when scaling.
    """
    if not shutil.which("signalp"):
        raise FileNotFoundError("signalp not found on system $PATH")

    fp = NTF("w", delete=False)
    fasta = make_fasta(proteins)
    output = None
    try:
        with fp:
            fp.write(fasta)
        params = {
            "path": "signalp",
            "-fasta": fp.name,
            "-batch": str(batch),
            "-format": "short",
            "last": "-stdout",
        }
        command = make_command(params)
        output = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            universal_newlines=True,
        ).stdout
    finally:
        os.unlink(fp.name)
    return output


def lfasta(sequence, matrix):
    """Run lfasta on protein sequence with given scoring matrix.

    RADAR combines output of BLOSUM50 and PAM250 matrices for
    repeat detection.

    This function must be called under a NTF()
    context manager.

    Args:
        sequence (str): path of tempfile containing
                                       amino acid sequence
        matrix (str): which scoring matrix to use
    Returns:
        (str): stdout from subprocess call
    """
    assert matrix in ["BLOSUM50", "PAM250"]

    command = [get_path("lfasta"), sequence, sequence]

    if matrix == "PAM250":
        command.append("-f -12 -g -2 -s 250")
    return subprocess.run(
        command, stdout=subprocess.PIPE, universal_newlines=True
    ).stdout


def radar(protein):
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

    with NTF() as blosum, NTF() as pam, NTF() as sequence:
        sequence.write(protein.sequence.encode())
        sequence.flush()

        # Run lfasta with BLOSUM50 matrix, then ensure we've gone back to start
        blosum.write(lfasta(sequence.name, "BLOSUM50").encode())
        sequence.seek(0)

        # Run lfasta with PAM250 matrix
        pam.write(lfasta(sequence.name, "PAM250").encode())

        blosum.seek(0)
        pam.seek(0)

        params = {
            "path": get_path("radar"),
            "-P": sequence.name,
            "-R": sequence.name,
            "-Q": blosum.name,
            "-S": pam.name,
        }

        return subprocess.run(
            make_command(params), stdout=subprocess.PIPE, universal_newlines=True
        ).stdout
