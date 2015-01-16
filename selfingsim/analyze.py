from __future__ import print_function

# standard imports
import argparse
from collections import Counter
import math

# within-package imports
from utils import getdata, getn, getnloc

def main():
    sp = argparse.ArgumentParser().add_subparsers()
    setup_command_line(sp)
    a = sp.parse_args()
    a.func(a)

def setup_command_line(sp):
    # setup shared command line argument
    pp = argarse.ArgumentParser(add_help = False)
    pp.add_argument(
            "samplefile",
            type = str,
            help = "file containing samples (output of (sub)sample command")
    pp.add_argument(
            "-n",
            acton = "store_true",
            default = False,
            help = "set this flag to suppress headder line")

    # setup command line arguments for indiviudal subcommands
    p = sp.add_parser("inbcoeff", parents = [pp])
    p.set_defaults(func = inbcoeff)

    p = sp.add_parser("inbtime", parents = [pp])
    p.set_defaults(func = inbtime)

def inbcoeff(a):
    data = getdata(a.samplefile)
    n = getn(data)

    nloc = getnloc(data)
    tdata = t(genes(data))
    hobs = gethobs(tdata)
    hexp = gethexp(tdata)
    nalleles = getnumalleles(tdata)
    printlong(not a.n, a.samplefile, hobs, hexp, nalleles, tdata, [n])

def inbtime(a):
    data = getdata(a.samplefile)

    if not a.n:
        print("\t".join(["dataset", "sample", "inbreeding.time"]))
    for d in data:
        print("{}\tsample.{}\t{}".format(a.samplefile, d[0], d[1]))

def eofratio(fis):
    data = [f for f in fis if not math.isnan(f)]
    return sum(data) / len(data) if len(data) > 0 else float('nan')

def ratioofes(hobs, hexp):
    ho = sum(hobs)
    he = sum(hexp)
    return 1.0 - ho / he if he != 0.0 else float('nan')

def ratioofes2(hobs, hexp, n):
    nn = 2 * n[0]
    ho = sum(hobs)
    he = sum(hexp)
    c = sum(h / nn for h in hobs)
    return (he - ho + c) / (he - c) if he - c != 0.0 else float('nan')


def printlong(hasheader, fname, hobs, hexp, nalleles, data, size):
    nn = 2 * size[0]
    if hasheader:
        print("\t".join(["dataset", "locus", "Hobs", "Hexp", "Fis", "Fis.corrected", "alleles"]))
    Fis = [1.0 - ho / he if he != 0.0 else float('nan') for ho, he in zip(hobs, hexp)]
    Fis2 = [(he - ho + ho / nn) / (he - ho / nn) if he - ho / nn != 0.0 else float('nan') for ho, he in zip(hobs, hexp)]
    for i, (j, k, F, F2, n) in enumerate(zip(hobs, hexp, Fis, Fis2, nalleles)):
        print("{}\tlocus.{}\t{}\t{}\t{}\t{}\t{}".format(fname, i, j, k, F, F2, n))
    print("{}\toverall.E.of.ratio\t{}\t{}\t{}\t{}\t{}".format(
        fname,
        sum(hobs) / len(hobs),
        sum(hexp) / len(hexp),
        eofratio(Fis),
        eofratio(Fis2),
        sum(nalleles) / float(len(nalleles))
        )
        )
    print("{}\toverall.ratio.of.Es\t{}\t{}\t{}\t{}\t{}".format(
        fname,
        sum(hobs) / len(hobs),
        sum(hexp) / len(hexp),
        ratioofes(hobs, hexp),
        ratioofes2(hobs, hexp, size),
        sum(nalleles) / float(len(nalleles))
        )
        )

    def genes(data):
        return [d[2] for d in data]

def t(data):
    return list(zip(*data))

def gethomofreq(data):
    size = len(data[0])
    count = [Counter(genot[0] for genot in locus if genot[0] == genot[1]).values() for locus in data]
    return [[float(val) / size for val in locus] for locus in count]

def getallelefreq(data):
    size = 2 * len(data[0])
    count = [Counter(allele for genot in locus for allele in genot).values() for locus in data]
    return [[float(val) / size for val in locus] for locus in count]

def gethobs(data):
    hs = gethomofreq(data)
    return [1.0 - sum(locus) for locus in hs]


def gethexp(data):
    als = getallelefreq(data)
    return [1.0 - sum(a * a for a in locus) for locus in als]

def getnumalleles(data):
    return [len(set(allele for genot in locus for allele in genot)) for locus in data]

if __name__ == '__main__':
    main()
