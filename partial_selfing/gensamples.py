from __future__ import print_function
import argparse, random, pickle, operator, json, sys, itertools, math, os.path
from collections import Counter, defaultdict

def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()
    ssp = sp.add_parser("sample")
    sssp = sp.add_parser("subsample")
    psp = sp.add_parser("phase")
    p2rsp = sp.add_parser("phase2rmes")
    nsp = sp.add_parser("nexus")
    rsp = sp.add_parser("rmes")
    asp = sp.add_parser("inbcoeff")
    isp = sp.add_parser("inbtime")
    rcsp = sp.add_parser("rmescombine")

    ssp.add_argument(
        "FILE",
        type = str
    )
    ssp.add_argument(
        "GEN",
        type = int
    )
    ssp.add_argument(
        "SAMPLE",
        type = int
    )
    ssp.add_argument(
        "REP",
        type = int
    )
    sssp.add_argument(
        "INPUT",
        type = str
    )
    sssp.add_argument(
        "OUTPUT",
        type = str
    )
    sssp.add_argument(
        "SAMPLE",
        type = int
    )
    sssp.add_argument(
        "LOCI",
        type = int,
        nargs = "*"
    )
    psp.add_argument(
        "FILE",
        type = str
    )
    p2rsp.add_argument(
        "FILE",
        type = str
    )
    p2rsp.add_argument(
        "-w",
        action = "store_true",
        default = False
    )
    nsp.add_argument(
        "FILE",
        type = str
    )
    nsp.add_argument(
        "SIZE",
        type = int,
        nargs = "*"
    )
    nsp.add_argument(
        "-w",
        action = "store_true",
        default = False
    )
    rsp.add_argument(
        "FILE",
        type = str
    )
    rsp.add_argument(
        "-w",
        action = "store_true",
        default = False
    )
    asp.add_argument(
        "FILE",
        type = str
    )
    asp.add_argument(
        "-n",
        action = "store_true",
        default = False
    )
    isp.add_argument(
        "FILE",
        type = str
    )
    isp.add_argument(
        "-n",
        action = "store_true",
        default = False
    )
    rcsp.add_argument(
        "OUT",
        type =str
    )
    rcsp.add_argument(
        "FILES",
        nargs = "+"
    )
    rcsp.add_argument(
        "-w",
        action = "store_true",
        default = False
    )

    ssp.set_defaults(func = sample)
    sssp.set_defaults(func = subsample)
    psp.set_defaults(func = phase)
    p2rsp.set_defaults(func = phase2rmes)
    nsp.set_defaults(func = nexus)
    rsp.set_defaults(func = rmes)
    asp.set_defaults(func = inbcoeff)
    isp.set_defaults(func = inbtime)
    rcsp.set_defaults(func = rmescombine)

    a = p.parse_args()
    a.func(a)

def sample(a):
    fbase = a.FILE.split(".")[:-1]
    data = getgeneration(a.FILE, a.GEN)

    for i in xrange(a.REP):
        s = getsample(data, a.SAMPLE)
        ofname = ".".join(
            fbase +
            [
                str(a.SAMPLE),
                "{:03}".format(i + 1),
                "json"
            ]
        )
        with open(ofname, "w") as wf:
            json.dump(s, wf)


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

def subsample(a):
    with open(a.INPUT, "r") as rf:
        data = json.load(rf)
        size = getn(data)
        loc = getnloc(data)
        if size < a.SAMPLE:
            print("SAMPLE too large", file = sys.stderr)
            sys.exit(1)
        if loc < len(a.LOCI):
            print("Too many loci", file = sys.stderr)

        samples = random.sample(data, a.SAMPLE)
        d = [[s[0], s[1], list(s[2][i] for i in a.LOCI)] for s in samples]
        with open(a.OUTPUT, "w") as wf:
            json.dump(d, wf)


def phase(a):
    ofname = ".".join(a.FILE.split(".")[:-1] + ["phase"])
    with open(a.FILE, "r") as rf, open(ofname, "w") as wf:
        data = json.load(rf)
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

    ofname = ".".join(a.FILE.split(".")[:-1] + ["nex"])
    with open(a.FILE, "r") as rf:
        data = json.load(rf)

    n = getn(data)
    if len(a.SIZE) > 0 and n != sum(a.SIZE):
        print("Sum of subsample sizes does not match the total sample size", file = sys.stderr)
        sys.exit(1)
    elif len(a.SIZE) > 0:
        size = a.SIZE
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
    print("Processing {}...".format(a.FILE))
    nl = getnlchar(a)

    ofname = ".".join(a.FILE.split(".")[:-1] + ["rmes"])
    with open(a.FILE, "r") as rf, open(ofname, "w") as wf:
        data = json.load(rf)
        size = getn(data)
        loc = getnloc(data)
        print(1, file = wf, end = nl)
        print(a.FILE, file = wf, end = nl)
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

def getnlchar(a):
    if a.w:
        return "\r\n"
    else:
        return "\n"



def inbcoeff(a):
    data = getdata(a.FILE)
    n = getn(data)

    nloc = getnloc(data)
    tdata = t(genes(data))
    hobs = gethobs(tdata)
    hexp = gethexp(tdata)
    nalleles = getnumalleles(tdata)
    printlong(not a.n, a.FILE, hobs, hexp, nalleles, tdata, [n])


def eofratio(fis):
    data = [f for f in fis if not math.isnan(f)]
    return sum(data) / len(data) if len(data) > 0 else float('nan')

def ratioofes(hobs, hexp):
    ho = sum(hobs)
    he = sum(hexp)
    return 1.0 - ho / he if he != 0.0 else float('nan')

def printlong(hasheader, fname, hobs, hexp, nalleles, data, size):
    if hasheader:
        print("\t".join(["dataset", "locus", "Hobs", "Hexp", "Fis", "alleles"]))
    Fis = [1.0 - ho / he if he != 0.0 else float('nan') for ho, he in zip(hobs, hexp)]
    for i, (j, k, F, n) in enumerate(zip(hobs, hexp, Fis, nalleles)):
        print("{}\tlocus.{}\t{}\t{}\t{}\t{}".format(fname, i, j, k, F, n))
    print("{}\toverall.E.of.ratio\t{}\t{}\t{}\t{}".format(
            fname,
            sum(hobs) / len(hobs),
            sum(hexp) / len(hexp),
            eofratio(Fis),
            sum(nalleles) / float(len(nalleles))
        )
    )
    print("{}\toverall.ratio.of.Es\t{}\t{}\t{}\t{}".format(
            fname,
            sum(hobs) / len(hobs),
            sum(hexp) / len(hexp),
            ratioofes(hobs, hexp),
            sum(nalleles) / float(len(nalleles))
        )
    )


def inbtime(a):
    data = getdata(a.FILE)

    if not a.n:
        print("\t".join(["dataset", "sample", "inbreeding.time"]))
    for d in data:
        print("{}\tsample.{}\t{}".format(a.FILE, d[0], d[1]))


def getdata(f):
    with open(f, "r") as rf:
        return json.load(rf)


def getn(data):
    return len(data)


def getnloc(data):
    if len(data) > 0:
        return len(data[0][2])
    else:
        return 0

def getsubboundaries(size):
    begin = [0]
    end = []
    for s in size[:-1]:
        begin.append(begin[-1] + s)
        end.append(begin[-1])
    end.append(sum(size))
    return begin, end

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

def rmescombine(a):
    nl = getnlchar(a)

    hds, bds = [], []

    npop = 0

    for f in a.FILES:
        with open(f, "r") as rf:
            next(rf)
            npop += 1
            for dummy in range(3):
                hds.append(next(rf).rstrip())
            for l in rf:
                bds.append(l.rstrip())

    with open(a.OUT, "w") as wf:
        print(npop, file = wf, end = nl)
        for h in hds:
            print(h, file = wf, end = nl)
        for b in bds:
            print(b, file = wf, end = nl)

def phase2rmes(a):
    popname = a.FILE
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

if __name__ == "__main__":
    main()
