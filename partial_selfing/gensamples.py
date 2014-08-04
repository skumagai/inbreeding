from __future__ import print_function
import argparse, random, pickle, operator, json

def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()
    ssp = sp.add_parser("sample")
    psp = sp.add_parser("phase")
    nsp = sp.add_parser("nexus")

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

    ssp.set_defaults(func = sample)
    psp.set_defaults(func = phase)
    nsp.set_defaults(func = nexus)

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
        size = len(data)
        loc = len(data[0][2])
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
    if a.w:
        nl = "\r\n"
    else:
        nl = "\n"

    ofname = ".".join(a.FILE.split(".")[:-1] + ["nex"])
    with open(a.FILE, "r") as rf, open(ofname, "w") as wf:
        data = json.load(rf)
        size = len(data)
        loc = len(data[0][2])
        print("#NEXUS", file = wf, end=nl)
        print("begin gdadata;", file = wf, end=nl)
        print("dimensions npops={} nloci={};".format(1, loc), file = wf, end=nl)
        print("format missing=? separator=/;", file = wf, end=nl)
        print("matrix", file = wf, end=nl)
        print("pop1:", file = wf, end=nl)
        for d in data:
            print(
                " ".join(
                    ["sample.{}".format(d[0])] +
                    ["{}/{}".format(*i) for i in d[2]]
                ),
                file = wf,
                end=nl
            )
        print(";", file = wf, end=nl)
        print("end;", file = wf, end=nl)


if __name__ == "__main__":
    main()
