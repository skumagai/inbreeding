from __future__ import print_function
import argparse, random, pickle, operator, json
from collections import Counter, defaultdict

def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()
    ssp = sp.add_parser("sample")
    psp = sp.add_parser("phase")
    nsp = sp.add_parser("nexus")
    rsp = sp.add_parser("rmes")
    asp = sp.add_parser("inbcoeff")
    isp = sp.add_parser("inbtime")
    rcsp = sp.add_parser("rmescombine")

    ssp.add_argument(
        "FILE",
        type = str,
    )
    ssp.add_argument(
        "GEN",
        type = int,
    )
    ssp.add_argument(
        "SAMPLE",
        type = int,
    )
    ssp.add_argument(
        "REP",
        type = int
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
    with open(a.FILE, "r") as rf, open(ofname, "w") as wf:
        data = json.load(rf)
        size = getn(data)
        loc = getnloc(data)
        print("#NEXUS", file = wf, end = nl)
        print("begin gdadata;", file = wf, end=nl)
        print("dimensions npops={} nloci={};".format(1, loc), file = wf, end = nl)
        print("format missing=? separator=/;", file = wf, end = nl)
        print("matrix", file = wf, end = nl)
        print("pop1:", file = wf, end = nl)
        for d in data:
            print(
                " ".join(
                    ["sample.{}".format(d[0])] +
                    ["{}/{}".format(*i) for i in d[2]]
                ),
                file = wf,
                end=nl
            )
        print(";", file = wf, end = nl)
        print("end;", file = wf, end = nl)


def rmes(a):
    nl = getnlchar(a)

    ofname = ".".join(a.FILE.split(".")[:-1] + ["rmes"])
    with open(a.FILE, "r") as rf, open(ofname, "w") as wf:
        data = json.load(rf)
        size = getn(data)
        loc = getnloc(data)
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
    hobs = gethobs(data, n, nloc)
    hexp = gethexp(data, n, nloc)
    f = getf(data, n, nloc)
    fall = getfall(data, n, nloc)
    printlong(not a.n, a.FILE, hobs, hexp, f, fall)


def printlong(hasheader, fname, hobs, hexp, f, fall):
    if hasheader:
        print("\t".join(["dataset", "locus", "Hobs", "Hexp", "f"]))

    for i, (j, k, l) in enumerate(zip(hobs, hexp, f)):
        print("{}\tlocus.{}\t{}\t{}".format(fname, i, j, k))
    print("{}\toverall\t{}\t{}".format(
            fname,
            sum(hobs) / len(hobs),
            sum(hexp) / len(hobs),
            fall
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
        return len(data[1][2])
    else:
        return 0

def gethobs(data, n, nloc):
    nf = 1.0 * n
    record = [0 for dummy in range(nloc)]
    for d in data:
        da = d[2]
        for i in range(nloc):
            record[i] += 1 if da[i][0] != da[i][1] else 0
    return [i / nf for i in record]


def gethexp(data, n, nloc):
    f = 2.0 * n / (2 * n - 1)
    ps = getp(data, n, nloc)
    return [f * (1.0 - sum([j**2.0 for j in i.values()])) for i in ps]


def getf(data, n, nloc):
    return [
        i / j if j != 0.0 else 0.0
        for i, j in zip(*getfnumdenom(data, n, nloc))
    ]


def getfall(data, n, nloc):
    num, denom = getfnumdenom(data, n, nloc)
    return sum(num) / sum(denom)


def getfnumdenom(data, n, nloc):
    nf2 = 2.0 * n
    ps = getp(data, n, nloc)
    Ps = getP(data, n, nloc)
    dPs = [defaultdict(lambda: 0.0, i) for i in Ps]

    tail = [
        (1.0 - sum([k for j, k in i.items() if j[0] == j[1]])) / nf2
        for i in dPs
    ]

    num = [
        sum(j[(l, l)] - m**2.0 for l, m in i.items()) + k
        for i, j, k in zip(ps, dPs, tail)
    ]

    denom = [
        (1.0 - sum(k**2.0 for k in i.values())) - j
        for i, j in zip(ps, tail)
    ]

    return num, denom


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

    for f in a.FILES:
        with open(f, "r") as rf:
            for dummy in range(3):
                hds.append(next(rf).rstrip())
            for l in rf:
                bds.append(l.rstrip())

    with open(a.OUT, "w") as wf:
        for h in hds:
            print(h, file = wf, end = nl)
        for b in bds:
            print(b, file = wf, end = nl)

if __name__ == "__main__":
    main()
