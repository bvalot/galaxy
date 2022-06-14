#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Return a filter fasta file for list of accession"""

import sys
import argparse
from Bio import SeqIO

desc = "Filter a fasta file based on ."
command = argparse.ArgumentParser(
    prog="fasta_filter_to_accession",
    description=desc,
    usage="%(prog)s [options] fasta accessions",
)
command.add_argument(
    "-o",
    "--output",
    default=sys.stdout,
    type=argparse.FileType("w"),
    nargs="?",
    help="Output result to, default:stdout",
)
command.add_argument("fasta", type=argparse.FileType("r"), help="Fasta file to filter")
command.add_argument(
    "access",
    metavar="accessions",
    type=str,
    nargs="+",
    help="List of accession to extract",
)

if __name__ == "__main__":
    """Performed job on execution script"""
    args = command.parse_args()
    output = args.output
    access = set(args.access)
    for seq in SeqIO.parse(args.fasta, "fasta"):
        if seq.id in access:
            SeqIO.write(seq, output, "fasta")
            access.remove(seq.id)
    if len(access) > 0:
        sys.stderr.write("Accessions not found : " + str(" - ").join(access) + "\n")
