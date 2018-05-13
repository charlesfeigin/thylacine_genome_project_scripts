#!/usr/bin/env python

from __future__ import print_function
import sys
import os, os.path
import numpy as np
import numpy.random as random

def parse_fai(filename, idir, suffix):
    scaffolds = []
    lengths = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            fields = line.split("\t")

            scaffold = fields[0]
            length = int(fields[1])

            # we deleted zero length files, so check for existence
            ifn = "{}/{}{}".format(idir, scaffold, suffix)
            if not os.path.exists(ifn):
                continue

            scaffolds.append(scaffold)
            lengths.append(length)

    return scaffolds, lengths

def parse_tsv(filename):
    with open(filename) as f:
        for line in f:
            line = line.strip()
            fields = line.split("\t")
            yield fields

def replace_col1(ifn, ofn, new_col1):
    with open(ofn, "w") as f:
        for fields in parse_tsv(ifn):
            ofields = [new_col1] + fields[1:]
            print(*ofields, sep="\t", file=f)

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("usage: {} in.fai idir odir".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    suffix = ".MSMC.txt"
    nreps = 100
    in_fai = sys.argv[1]
    idir = sys.argv[2]
    odir = sys.argv[3]

    scaffolds, lengths = parse_fai(in_fai, idir, suffix)
    total_length = np.sum(lengths)

    for i in range(nreps):
        opath = "{}/bootstrap_{}".format(odir, i)
        if not os.path.exists(opath):
            os.makedirs(opath)
        length = 0
        n = 0
        while length < total_length:
            j = random.randint(len(scaffolds))
            scaffold = scaffolds[j]
            length += lengths[j]
            n += 1
            ifn = "{}/{}{}".format(idir, scaffold, suffix)
            ofn = "{}/{}{}".format(opath, n, suffix)
            replace_col1(ifn, ofn, str(n))
