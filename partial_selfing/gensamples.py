from __future__ import print_function
import argparse, random, pickle, operator, json, sys, itertools
from collections import Counter, defaultdict

def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()
    ssp = sp.add_parser("sample")
    sssp = sp.add_parser("subsample")
    psp = sp.add_parser("phase")
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
    pdata = partitiondata(data, [n])
    hobs = gethobs(pdata, [n])
    hexp = gethexp(pdata, [n])
    printlong(not a.n, a.FILE, hobs, hexp, pdata, [n])


def printlong(hasheader, fname, hobs, hexp, data, size):
    if hasheader:
        print("\t".join(["dataset", "locus", "Hobs", "Hexp", "Fis", "f", "g"]))
    Fis = [1.0 - ho / he for ho, he in zip(hobs, hexp)]
    for i, (j, k, F) in enumerate(zip(hobs, hexp, Fis)):
        print("{}\tlocus.{}\t{}\t{}\t{}".format(fname, i, j, k, F))
    print("{}\toverall\t{}\t{}\t{}".format(
            fname,
            sum(hobs) / len(hobs),
            sum(hexp) / len(hobs),
            sum(Fis) / len(Fis),
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

def transposedata(data):
    return [list(v) for v in zip(*[ind[2] for ind in data])]

def partitiondata(data, size):
    begin, end = getsubboundaries(size)
    nal = [float(d) for d in size]
    # tranpose data: now row (outer-most list) corresponds to locus rather than individual
    tdata = transposedata(data)
    # partition data into per-subpopulation basis
    return [[list(itertools.islice(locus, b, e)) for b, e in zip(begin, end)] for locus in tdata]

def gethomofreq(data, size):
    count = [[Counter(genot[0] for genot in deme if genot[0] == genot[1]).values() for deme in locus] for locus in data]
    return [[[float(val) / s for val in deme] for deme, s in zip(locus, size)] for locus in count]

def getallelefreq(data, size):
    count = [[Counter(allele for genot in deme for allele in genot).values() for deme in locus] for locus in data]
    return [[[float(val) / (2 * s) for val in deme] for deme, s in zip(locus, size)] for locus in count]

def gethobs(data, size):
    hs = gethomofreq(data, size)
    freqs = [float(s) / sum(size) for s in size]
    return [1.0 - sum(sum(deme) * f for deme, f in zip(locus, freqs)) for locus in hs]


def gethexp(data, size):
    als = getallelefreq(data, size)
    freqs = [float(s) / sum(size) for s in size]
    return [1.0 - sum(sum(a * a for a in deme) * f for deme, f in zip(locus, freqs)) for locus in als]


def unstructuredata(data):
    return [list(itertools.chain(*locus)) for locus in data]

def getf(data, size):
    fs, gs = getfg(data, size)
    return [(f - g) / (1.0 - g) for f, g in zip(fs, gs)]


def getfg(data, size):
    nf = float(sum(size))
    data = unstructuredata(data)
    fs = [len(list(True for ind in locus if ind[0] == ind[1])) / nf for locus in data]
    gdata = unstructuredata(data)
    ng = len(gdata[1])
    gs = [float(len(list(True for i in range(ng) for j in range(i + 1, ng) if locus[i] == locus[j]))) / (ng * (ng - 1) / 2) for locus in gdata]
    return fs, gs

def getfall(data, size):
    f = getf(data, size)
    return sum(f) / len(f)

def getp(data, n, nloc):
    nf2 = 2 * float(n)
    a = [[] for i in range(nloc)]
    for d in data:
        da = d[2]
        for i in range(nloc):
            a[i].extend(da[i])
    return [{j: k / nf2 for j, k in Counter(i).items()} for i in a]


def getP(data, n, nloc):
    nf = float(n)
    g = [[] for i in range(nloc)]
    for d in data:
        da = d[2]
        for i in range(nloc):
            g[i].append(tuple(da[i]))
    return [{j: k / nf for j, k in Counter(i).items()} for i in g]


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

if __name__ == "__main__":
    main()
