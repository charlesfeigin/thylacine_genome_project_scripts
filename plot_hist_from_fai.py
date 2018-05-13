#!/usr/bin/env python

from __future__ import print_function
import sys
import os.path
import matplotlib
matplotlib.use('Agg') # don't try to use $DISPLAY
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.gridspec as gridspec
import numpy as np

def parse_fai(filename):
    l = []
    with open(filename) as f:
        for line in f:
            line = line.strip()
            fields = line.split("\t")

            size = int(fields[1])
            l.append(size)

    return l


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: {} in.fai out.pdf".format(sys.argv[0]), file=sys.stderr)
        exit(1)

    l = np.array(parse_fai(sys.argv[1]), dtype=float)
    l /= 1e6

    title = os.path.basename(sys.argv[1])
    if title.endswith(".fai"):
        title = title[:-4]

    pdf = PdfPages(sys.argv[2])
    #fig_w, fig_h = plt.figaspect(9.0/16.0)
    fig_w, fig_h = plt.figaspect(3.0/4.0)
    fig1 = plt.figure(figsize=(fig_w, fig_h))
    gs1 = gridspec.GridSpec(1, 1)
    ax1 = fig1.add_subplot(gs1[0])

    ax1.hist(l, bins=30)
    ax1.set_xlabel("scaffold length (mb)")
    ax1.set_ylabel("counts")
    ax1.set_title(title)


    plt.tight_layout()
    pdf.savefig()
    pdf.close()
