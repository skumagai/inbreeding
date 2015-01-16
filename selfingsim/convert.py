# standard imports
import argparse
import itertools
import os.path

# within-package imports
from utils import getdata, getn, getnloc

def main():
    sp = argparse.ArgumentParser().add_subparsers()
    setup_commnad_line(sp)
    a = sp.parse_args()
    a.func(a)

def setup_command_line(sp):
    # setup shared command line arguments
    pp1 = argparse.ArgumentParser(add_help = False)
    pp1.add_argument(
            "samplefile",
            type = str,
            help = "file containing samples (output of (sub)sample command)")

    pp2 = argparse.ArgumentParser()
    pp2.add_argument(
            "-w",
            action = "store_true",
            default = False,
            help = "set this flag to force Windows newline character (CRLF).")

    # setup command line arguments for individual subcommands
    p = sp.add_parser("phase", parents = [pp1])
    p.set_defaults(func = phase)

    p = sp.add_parser("phase2rmes", parents = [pp2])
    p.add_argument(
            "phasefile",
            type = str,
            help = "file containing samples in phase format.")
    p.set_defaults(func = phase2rmes)

    p = sp.add_parser("nexus", parents = [pp1, pp2])
    p.add_argument(
            "localsizes",
            type = int,
            nargs = "*",
            help = "local sample sizes (zero or more)"
            )
    p.set_defaults(func = nexus)

    p = sp.add_parser("rmes", parents = [pp1, pp2])
    p.set_defaults(func = rmes)

    p = sp.add_parser("rmescombine", parents = [pp2])
    p.add_argument(
            "combinedfile",
            type =str,
            help = "file of output")
    p.add_argument(
            "rmesfiles",
            type = str,
            nargs = "+",
            help = "one or more sample files in RMES format")
    p.set_defaults(func = rmescombine)

def phase(a):
    ofname = ".".join(a.samplefile.split(".")[:-1] + ["phase"])
    data = getdata(a.samplefile)
    with open(ofname, "w") as wf:
        size = getn(data)
        loc = getnloc(data)
        print(
                "{}\n{}\n{}".format(
                    size,
                    loc,
                    "M" * loc
                    ),
                file = wf
                )
        for d in data:
            print(
                    "sample.{}\t{}".format(
                        d[0],
                        "\t".join([str(j) for i in d[2] for j in i])
                        ),
                    file = wf
                    )

            def nexus(a):
                nl = getnlchar(a)

    ofname = ".".join(a.samplefile.split(".")[:-1] + ["nex"])
    data = getdata(a.samplefile)

    n = getn(data)
    if len(a.localsizes) > 0:
        if n != sum(a.localsizes):
            print("Sum of subsample sizes does not match the total sample size", file = sys.stderr)
            sys.exit(1)
        size = a.localsizes
    else:
        size = [n]

    with open(ofname, "w") as wf:
        loc = getnloc(data)
        print("#NEXUS", file = wf, end = nl)
        print("begin gdadata;", file = wf, end=nl)
        print("dimensions npops={} nloci={};".format(len(size), loc), file = wf, end = nl)
        print("format missing=? separator=/;", file = wf, end = nl)
        print("matrix", file = wf, end = nl)
        begin, end = getsubboundaries(size)
        subs = [itertools.islice(data, b, e) for b, e in zip(begin, end)]
        for i, sub in enumerate(subs):
            print("subpop.{}:".format(i), file = wf, end = nl)
            for d in sub:
                print(
                        " ".join(
                            ["sample.{}".format(d[0])] +
                            ["{}/{}".format(*j) for j in d[2]]
                            ),
                        file = wf,
                        end=nl
                        )
                if i < len(size) - 1:
                    print(",", file = wf, end = nl)
        print(";", file = wf, end = nl)
        print("end;", file = wf, end = nl)

def rmes(a):
    print("Processing {}...".format(a.samplefile))
    nl = getnlchar(a)

    ofname = ".".join(a.samplefile.split(".")[:-1] + ["rmes"])
    data = getdata(a.samplefile)
    with open(ofname, "w") as wf:
        size = getn(data)
        loc = getnloc(data)
        print(1, file = wf, end = nl)
        print(a.samplefile, file = wf, end = nl)
        print(size, file = wf, end = nl)
        print(loc, file = wf, end = nl)
        for d in data:
            print(
                    " ".join([
                        "0" if g[0] == g[1] else "1"
                        for g in d[2]
                        ]),
                    file = wf,
                    end = nl
                    )

            def rmescombine(a):
                nl = getnlchar(a)

    hds, bds = [], []

    npop = 0

    for f in a.rmesfiles:
        with open(f, "r") as rf:
            next(rf)
            npop += 1
            for dummy in range(3):
                hds.append(next(rf).rstrip())
            for l in rf:
                bds.append(l.rstrip())

    with open(a.combinedfile, "w") as wf:
        print(npop, file = wf, end = nl)
        for h in hds:
            print(h, file = wf, end = nl)
        for b in bds:
            print(b, file = wf, end = nl)

def phase2rmes(a):
    popname = a.phasefile
    outfile = os.path.splitext(popname)[0] + ".rmes"
    nl = getnlchar(a)
    data = []
    with open(popname, "r") as rf:
        norg = next(rf).strip()
        nloc = next(rf).strip()
        next(rf)
        for l in rf:
            it = iter(l.strip().split("\t")[1:])
            cols = []
            for i, j in itertools.izip(it, it):
                if i == "NA" or j == "NA":
                    cols.append("-9")
                elif i == j:
                    cols.append("0")
                else:
                    cols.append("1")

            data.append(" ".join(cols))

    with open(outfile, "w") as wf:
        print(1, file = wf, end = nl)
        print(popname, file = wf, end = nl)
        print(norg, file = wf, end = nl)
        print(nloc, file = wf, end = nl)
        for d in data:
            print(d, file = wf, end = nl)

def getnlchar(a):
    if a.w:
        return "\r\n"
    else:
        return "\n"

def getsubboundaries(size):
    begin = [0]
    end = []
    for s in size[:-1]:
        begin.append(begin[-1] + s)
        end.append(begin[-1])
    end.append(sum(size))
    return begin, end
