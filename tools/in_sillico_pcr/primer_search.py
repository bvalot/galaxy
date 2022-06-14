#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""Search primers on database and return result"""

import argparse
import os
import subprocess
import sys
import tempfile

from Bio import SeqIO

gasst_exe = "Gassst"
nucleo = {
    "R": ("A", "G"),
    "Y": ("C", "T"),
    "S": ("G", "C"),
    "W": ("A", "T"),
    "K": ("G", "T"),
    "M": ("A", "C"),
    "B": ("C", "G", "T"),
    "D": ("A", "G", "T"),
    "H": ("A", "C", "T"),
    "V": ("A", "C", "G"),
    "N": ("A", "T", "C", "G"),
}

desc = "Search primer on database and return amplicons."
command = argparse.ArgumentParser(
    prog="primer_search.py",
    description=desc,
    usage="%(prog)s [options] forward reverse database",
)
command.add_argument(
    "-o",
    "--out",
    nargs="?",
    type=argparse.FileType("w"),
    default=sys.stdout,
    help="Return amplicon on fasta file, default=stdout",
)
command.add_argument(
    "-d",
    "--debug",
    nargs="?",
    type=argparse.FileType("w+"),
    default=None,
    help="Conserve a copy of primer search, default No",
)
command.add_argument(
    "-e",
    "--error",
    nargs="?",
    type=int,
    default=1,
    help="Maximun error allowed on match, default=1",
)
command.add_argument(
    "-m",
    "--min",
    nargs="?",
    type=int,
    default=100,
    help="Min len amplicon size, default=100",
)
command.add_argument(
    "-M",
    "--max",
    nargs="?",
    type=int,
    default=1500,
    help="Max len amplicon size, default=1500",
)
command.add_argument(
    "-k",
    "--keep",
    action="store_true",
    help="Keep description instead of report PCR position",
)
command.add_argument(
    "-r", "--remove", action="store_true", help="Remove primer from reported amplicons"
)
command.add_argument("forward", type=str, help="Forward primer sequence")
command.add_argument("reverse", type=str, help="Reverse primer sequence")
command.add_argument(
    "database", type=argparse.FileType("r"), help="Database to search on fasta"
)
command.add_argument("-v", "--version", action="version", version="%(prog)s 0.3.0")


class GassstResult:
    """A simple Gassst class containing result
    rev : True indicate if reverse primers, else forward
    """

    def __init__(self, line):
        h = line.strip().split("\t")
        if len(h) != 9:
            raise Exception("Line seems not correspond to gassst result\n" + line)
        self.primer = h[0]
        self.rev = h[0] != "Forward"
        self.match = h[1]
        self.strand = int(h[2])
        self.start = int(h[3])
        self.error = int(h[5])
        self.stop = self.start + len(h[7])

    def __repr__(self):
        return (
            self.match
            + " : "
            + str(self.start)
            + " in strand:"
            + str(self.strand)
            + " in primer: "
            + self.primer
            + " with error: "
            + str(self.error)
        )

    def __eq__(self, other):
        return (
            isinstance(other, GassstResult)
            and self.primer == other.primer
            and self.rev == other.rev
            and self.match == other.match
            and self.strand == other.strand
            and self.start == other.start
        )

    def __ne__(self, other):
        return not self == other

    def __lt__(self, other):
        if self != other:
            raise Exception("Only compare equal GassstResult")
        return self.error < other.error


class Amplicon:
    """A simple Amplicon class containing forward and reverse gassst result"""

    def __init__(self, forw, reve):
        if forw.strand == reve.strand:
            raise Exception(
                "Amplicon could be in inverse strand for reverse/forwar primer\n"
                + forw
                + reve
            )
        self.forw = forw
        self.reve = reve
        self.strand = forw.strand

    def get_start(self):
        if self.strand == 0:
            return self.forw.start
        return self.reve.start

    def get_stop(self):
        if self.strand == 1:
            return self.forw.stop
        return self.reve.stop

    def __len__(self):
        return self.get_stop() - self.get_start()

    def __repr__(self):
        return (
            str(self.get_start())
            + " : "
            + str(self.get_stop())
            + " in strand:"
            + str(self.strand)
            + "; error F:"
            + str(self.forw.error)
            + ",R:"
            + str(self.reve.error)
        )


def degenerate_primer(seq):
    """Iterate all possible sequence
    from degenerate primer
    """
    notdegene = True
    for i, n in enumerate(seq.upper()):
        if n in nucleo:
            notdegene = False
            for mod in nucleo.get(n):
                seq_cut = list(seq)
                seq_cut[i] = mod
                for seq2 in degenerate_primer("".join(seq_cut)):
                    yield (seq2)
            break
    if notdegene:
        yield (seq)


def write_primer_fasta(forward_de, reverse_de):
    primer_fasta = tempfile.NamedTemporaryFile(delete=False, mode="w+t")
    for forward in degenerate_primer(forward_de):
        primer_fasta.write(">Forward\n" + forward + "\n")
    for reverse in degenerate_primer(reverse_de):
        primer_fasta.write(">Reverse\n" + reverse + "\n")
    primer_fasta.close()
    return primer_fasta.name


def performed_gassst(command):
    proc = subprocess.Popen(command, stderr=subprocess.PIPE, stdout=sys.stderr)
    error = ""
    for line in iter(proc.stderr.readline, b""):
        error += line.decode()
    if error != "":
        sys.stdout.write("Error during processing Gassst\n")
        raise Exception(error)


def read_result(f_res, f_debug):
    result = {}
    for line in iter(f_res.readline, ""):
        if line[0] == "F" or line[0] == "R":
            if f_debug is not None:
                f_debug.write(line)
            res = GassstResult(line)
            other_res = result.setdefault(res.match, {}).setdefault(res.primer, [])
            notfound = True
            for i, res2 in enumerate(other_res):
                if res == res2:
                    notfound = False
                    if res < res2:
                        other_res[i] = res
                        # print res
                        break
            if notfound:
                other_res.append(res)
    return result


def get_amplicons(result):
    """Function to get amplicon on a sequence from forward / reverse result"""
    for forw in result["Forward"]:
        for reve in result["Reverse"]:
            if forw.strand != reve.strand:
                amp = Amplicon(forw, reve)
                if amp.get_start() < amp.get_stop():  # must be side by side
                    yield amp
                else:
                    del amp


if __name__ == "__main__":
    """Performed job on execution script"""
    args = command.parse_args()
    primer_fasta = write_primer_fasta(args.forward, args.reverse)
    output_result = tempfile.NamedTemporaryFile(delete=False, mode="w+")
    identity = int(
        100 - args.error / float(max(len(args.forward), len(args.reverse))) * 100
    )

    command = [
        gasst_exe,
        "-p",
        str(identity),
        "-g",
        "0",
        "-w",
        "6",
        "-h",
        "0",
        "-i",
        primer_fasta,
        "-d",
        args.database.name,
        "-o",
        output_result.name,
    ]
    performed_gassst(command)
    result = read_result(output_result, args.debug)

    # load amplicons
    sys.stderr.write("Search amplicons : ")
    amplicons = {}
    for seq in result.keys():
        if len(result[seq]) == 2:
            # print "Search amplicon on seq: " + seq
            for amplicon in get_amplicons(result[seq]):
                if len(amplicon) > args.min and len(amplicon) < args.max:
                    amplicons.setdefault(seq, []).append(amplicon)
                else:
                    del amplicon
            # print "Find amplicons : "  + str(len(amplicons[seq]))
    sys.stderr.write(str(sum(map(len, amplicons.values()))) + "\n\n")

    # print amplicons

    sys.stderr.write("Write result : ")
    for record in SeqIO.parse(args.database, "fasta"):
        if record.id not in amplicons:
            continue
        for i, amplicon in enumerate(amplicons.get(record.id)):
            sub = record[amplicon.get_start() - 1: amplicon.get_stop() - 1]
            if amplicon.strand == 1:
                sub = sub.reverse_complement()
            if args.remove:
                sub = sub[len(args.forward): len(sub) - len(args.reverse)]
            if args.keep:
                sub.id = record.id
                sub.description = record.description
            else:
                sub.id = record.id + "_ampli_" + str(i)
                sub.description = str(amplicon)
            SeqIO.write(sub, args.out, "fasta")

    sys.stderr.write("OK\n")
    # delete tempfile
    os.unlink(primer_fasta)
    os.unlink(output_result.name)
