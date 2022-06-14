#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Return a filter fasta file for region of accession"""

import sys
import argparse
from Bio import SeqIO

desc = "Filter a fasta file based on region."
command = argparse.ArgumentParser(
    prog="fasta_filter_to_region",
    description=desc,
    usage="%(prog)s [options] fasta accessions:start-stop",
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
    help="List of region to extract, ex: accession:start-stop",
)

if __name__ == "__main__":
    """Performed job on execution script"""
    args = command.parse_args()
    output = args.output
    access = {}
    for acc in args.access:
        access.setdefault(acc.split(":")[0], []).append(acc.split(":")[1].split("-"))

    for seq in SeqIO.parse(args.fasta, "fasta"):
        accs = access.get(seq.id)
        if seq.id not in access:
            continue
        for acc in access.get(seq.id):
            sub = None
            if int(acc[0]) > int(acc[1]):
                # reverse
                sub = seq[int(acc[1]) - 1: int(acc[0])]
                sub.seq = sub.seq.reverse_complement()
                sub.id = seq.id + "_" + acc[0] + "_" + acc[1] + "_reverse"
            else:
                sub = seq[int(acc[0]) - 1: int(acc[1])]
                sub.id = seq.id + "_" + acc[0] + "_" + acc[1]
            SeqIO.write(sub, output, "fasta")
        access.pop(seq.id)
    if len(access) > 0:
        sys.stderr.write("Accessions not found : " + str(" - ").join(access) + "\n")
