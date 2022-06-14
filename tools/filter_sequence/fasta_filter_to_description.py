#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Return a filter fasta file containing key words on description"""

import argparse
import sys

from Bio import SeqIO

desc = "Filter a fasta file based on description."
command = argparse.ArgumentParser(
    prog="fasta_filter_to_description", description=desc, usage="%(prog)s [options] key"
)
command.add_argument(
    "-i",
    "--input",
    default=sys.stdin,
    type=argparse.FileType("r"),
    nargs="?",
    help="Input fasta file to, default:stdin",
)
command.add_argument(
    "-o",
    "--output",
    default=sys.stdout,
    type=argparse.FileType("w"),
    nargs="?",
    help="Output result to, default:stdout",
)
# command.add_argument('-e','--perl', action='store_true', \
#     help='Key word is an perl expression')
command.add_argument(
    "-V",
    "--reverse",
    action="store_true",
    help="Return sequence without key in description",
)
command.add_argument(
    "key", type=str, help="Key words that must be present in description"
)

if __name__ == "__main__":
    """Performed job on execution script"""
    args = command.parse_args()
    output = args.output
    count = 0
    for seq in SeqIO.parse(args.input, "fasta"):
        valid = False
        if args.reverse:
            if args.key not in seq.description:
                valid = True
        else:
            if args.key in seq.description:
                valid = True
        if valid:
            SeqIO.write(seq, output, "fasta")
            count += 1
    sys.stderr.write("Number of filter sequences : " + str(count) + "\n")
