from __future__ import print_function

# standard imports
import argparse
import json
import random
import sys

def main():
    "Run this script as stand-alone."
    p = argparse.ArgumentParser()
    setup_command_line(sp)
    args = p.parse_args()
    args.func(args)

def setup_command_line(sp):
    p = sp.add_parser("sample")
    p = sp.add_parser("subsample")

    p.add_argument(
            "simfile",
            type = str,
            help = "file containing simulation results")
    p.add_argument(
            "generations",
            type = int,
            help = "sampling generation")
    p.add_argument(
            "samplesize",
            type = int,
            help = "sample size")
    p.add_argument(
            "reps",
            type = int,
            "number of replicates")
    p.set_defaults(func = sample)

    p.add_argument(
            "samplefile",
            type = str,
            help = "file containing sample (output of sample command)")
    p.add_argument(
            "subsamplefile",
            type = str,
            help = "file storing subsamples")
    p.add_argument(
            "samplesize",
            type = int,
            help = "sample size")
    p.add_argument(
            "sampleloci",
            type = int,
            nargs = "*",
            help = "indicies of loci to sample (0-based index)")
    p.set_defaults(func = subsample)

def sample(a):
    fbase = a.simfile.split(".")[:-1]
    data = getgeneration(a.simfile, a.generations)
    digits = 1
    reps = float(a.reps)
    while reps > 10.:
        reps /= 10.
        digits += 1

    for i in xrange(a.reps):
        s = getsample(data, a.samplesize)
        ofname = ".".join(
                fbase +
                [
                    str(a.samplesize),
                    ("{:0" + str(digits) + "}").format(i + 1),
                    "json"
                    ]
                )
        with open(ofname, "w") as wf:
            json.dump(s, wf)

def subsample(a):
    with open(a.samplefile, "r") as rf:
        data = json.load(rf)
        size = getn(data)
        loc = getnloc(data)
        if size < a.samplesize:
            print("sub sample size too large", file = sys.stderr)
            sys.exit(1)
        if loc < len(a.sampleloci):
            print("Too many loci", file = sys.stderr)

        samples = random.sample(data, a.samplesize)
        d = [[s[0], s[1], list(s[2][i] for i in a.sampleloci)] for s in samples]
        with open(a.subsamplefile, "w") as wf:
            json.dump(d, wf)

def getgeneration(fname, gen):
    data = {}
    with open(fname, "r") as f:
        next(f)
        for line in f:
            fields = [int(round(float(i))) for i in line.split("\t")]
            if fields[1] == gen:
                try:
                    data[fields[0]]
                except KeyError:
                    data[fields[0]] = {}
                try:
                    data[fields[0]][fields[2]]["geno"].append(fields[5:])
                except KeyError:
                    data[fields[0]][fields[2]] = {
                            "selfing": fields[3],
                            "geno": [fields[5:]]
                            }
                    return data


def getsample(d, s, repidx = 0):
    size = len(d[repidx])
    return sorted(
            [simplify(i, d[repidx][i]) for i in random.sample(xrange(size), s)],
            key = operator.itemgetter(0)
            )

    def simplify(i, vals):
        return [
                i,
                vals["selfing"],
                list(
                    (i, j) for i, j in zip(*vals["geno"])
                    )
                ]

        if __name__ == '__main__':
            main()
