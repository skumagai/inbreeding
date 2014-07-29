import argparse, random, pickle, operator, json

def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()
    ssp = sp.add_parser("sample")
    psp = sp.add_parser("phase")

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

    ssp.set_defaults(func = sample)
    psp.set_defaults(func = phase)

    a = p.parse_args()
    a.func(a)

def sample(a):
    fbase = a.FILE.split(".")[:-1]
    data = getgeneration(a.FILE, a.GEN)

    for i in xrange(a.REP):
        s = getsample(data, a.SAMPLE)
        with open(".".join(fbase + ["{:03}".format(i), "json"]), "w") as wf:
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
        wf.write(
            "{}\n{}\n{}\n".format(
            size,
            loc,
            "M" * loc
        ))
        for d in data:
            wf.write("sample.{}\t{}\n".format(
                d[0],
                "\t".join([str(j) for i in d[2] for j in i])
            ))


if __name__ == "__main__":
    main()
