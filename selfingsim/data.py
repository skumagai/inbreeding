# standeard imports
from collections import Counter
import csv
from itertools import groupby, izip
import json
import math
import os.path
import random

def createsample(fname, gen = None):
    """
    Depending on file suffix, create an instance of an appropriate Sample object.

    The second argument "gen" is only meaningful when an instance is created from
    tsv file.
    """
    suffix = os.path.splitext(fname)[1]
    if suffix == ".tsv": # original simulation result
        return FullSample.fromtsv(fname, gen)
    elif suffix == ".json":
        return FullSample.fromjson(fname)
    elif suffix == ".phase":
        return BasicSample.fromphase(fname)
    raise ValueError


class BasicSample(object):
    """
    This class holds genotypic information of sample, and it also provides
    serialization to various formats (currently RMES, PHASE, and NEXUS).
    """
    def __init__(self, src, ids, genos):
        self._src = src
        self._ids = ids
        self._genos = genos
        self._nsam = len(genos)
        self._nloc = len(genos[0])

    @property
    def source(self):
        return self._src

    @property
    def ids(self):
        return self._ids

    @property
    def genotypes(self):
        return self._genos

    @property
    def n(self):
        return self._nsam

    @property
    def nloc(self):
        return self._nloc

    def sample(self, nsam, locs = None):
        """
        Sample a subset of individuals from the current BasicSample object.

        Optionally, you can specify which loci to include in the new sample.
        """
        locs = self._checklocs(locs)
        idx = self._drawindividuals(nsam)

        ids = [self._ids[i] for i in idx]
        genos = [[self._genos[i][j] for j in locs] for i in idx]

        return BasicSample(self._src, ids, genos)

    def _checklocs(self, locs):
        """
        Perform a basic sanity check if specified loci are a subset of existing
        loci. Additionally, it assumes all loci if no locus is specified.

        If this condition holds, this function returns a list of specified loci.
        Else raise ValueError.
        """
        # Sanity check: test if specified loci all exist.
        if locs == None:
            ret = range(len(self._genos[0]))
        else:
            nloc = len(self._geno[0])
            if len([i >= nloc for i in locs]) > 0:
                raise ValueError("Specified locus outside of range")
            ret = locs
        return ret


    def _drawindividuals(self, n):
        """
        Draw n individuals without replacement from the current sample.
        """
        return random.sample(xrange(self._nsam), n)

    @staticmethod
    def fromphase(fname):
        """
        Create an instance of BasicSample.
        """
        with open(fname, "rU") as f:
            nsam = int(next(f).strip())
            nloc = int(next(f).strip())
            types = next(f).strip()

            # Sanity check
            if len(types) != nloc:
                raise ValueError("Inconsistent number of loci.")

            ids = []
            genos = []

            for line in f:
                elems = line.strip().split()

                ids.append(elems[0])

                ielems = iter(elems[1:])
                genos.append([[i, j] for i, j in izip(ielems, ielems)])

        return BasicSample(fname, ids, genos)

    def tophase(self):
        """
        Write this sample as a phase-formatted string.
        """

        lines = [str(self._nsam), str(self._nloc), "M" * self._nloc]

        for idx, geno in zip(self._ids, self._genos):
            line = "\t".join([str(idx)] + [str(j) for i in geno for j in i])
            lines.append(line)

        return lines

    def tonexus(self, miss="?", sep = "/"):
        """
        Write this sample as a nexus-formatted string.
        """
        return tonexus([self], miss, sep)

    def tormes(self):
        """
        Write this sample as a RMES-formatted string.
        """
        return tormes([self])

    def inbreedingcoefficient(self):
        """
        Compute per-locus and overall inbreeding coefficient Fis.

        In addition to Fis, this method also returns other related statistics
        (observed and expected heterozygosities, and bias-corrected Fis as well as
        number of observed alleles.
        """
        genos = [[ind[loc] for ind in self._genos]
                for loc in xrange(self._nloc)]
        total = float(self.n)
        hobs = [sum(1 if g[0] != g[1] else 0 for g in geno) / total for geno in genos]

        alleles = [[gene for ind in self._genos for gene in ind[loc]]
                for loc in xrange(self._nloc)]
        total = float(2 * self.n)
        afreq = [[v / total for v in Counter(a).values()] for a in alleles]
        hexp = [1. - sum(a * a for a in locus) for locus in afreq]
        na = [len(a) for a in afreq]

        fis = [1. - ho / he if he != 0. else float('nan')
            for ho, he in zip(hobs, hexp)]
        fis2 = [(he - ho + ho / total) / (he - ho / total)
            if he - ho / total != 0.0 else float('nan')
            for ho, he in zip(hobs, hexp)]

        ho = sum(hobs)
        he = sum(hexp)
        c = ho / total

        fis.append((he - ho) / he if he != 0. else float('nan'))
        fis2.append((he - ho + c) / (he - c) if he - c != 0. else float('nan'))

        hobs.append(sum(hobs) / float(len(hobs)))
        hexp.append(sum(hexp) / float(len(hexp)))
        na.append(sum(na) / float(len(na)))

        keys = [i for i in xrange(self.nloc)]
        keys.append("overall")

        data = []
        for k, ho, he, n, f, f2 in zip(keys, hobs, hexp, na, fis, fis2):
            data.append({
                "key": k,
                "hobs": ho,
                "hexp": he,
                "fis": f,
                "fisc": f2,
                "nalleles": n})

        return data

class FullSample(BasicSample):
    """
    In addition to what BasicSample holds, this class also holds the number
    of generations until the last outcrossing event.
    """
    def __init__(self, src, ids, genos, inbgens):
        super(FullSample, self).__init__(src, ids, genos)
        self._inbgens = inbgens

    @property
    def tselfing(self):
        return self._inbgens

    def sample(self, nsam, locs = None):
        """
        Sample subset of individuals from the current sample.  Optionally,
        you can specify subsets of loci.
        """
        locs = self._checklocs(locs)
        idx = self._drawindividuals(nsam)

        ids = [self._ids[i] for i in idx]
        inbs = [self._inbgens[i] for i in idx]
        genos = [[self._genos[i][j] for j in locs] for i in idx]

        return FullSample(self._src, ids, genos, inbs)

    @staticmethod
    def fromtsv(fname, gen):
        """
        Create list of instances of FullSample from simulation results.

        Because instances are directly created from simulation results, they
        represent the entire population at the specified generation.

        If one file contains more than one replicates, this function returns
        separate FullSample instance for each of replicates.  If it records a
        single result, this function returns a list of length 1.

        If specified generations are not recorded, this returns an empty list.
        """
        samples = []
        with open(fname, "r") as f:
            reader = csv.reader(f, delimiter = "\t")
            # Throw out a header row
            next(f)
            rows = [row for row in reader if int(row[1]) == gen]

            # sanity check 1: total number of chromosomes (rows) must be multiple
            # of 2, as diploids have two chromosomes.
            if len(rows) % 2 != 0:
                raise ValueError("Number of chromosomes not mulitple of 2")


            for rep, rawdata in groupby(rows, lambda x: x[0]):
                irawdata = iter(rawdata)

                ids = []
                inbgens = []
                genos = []
                for v1, v2 in izip(irawdata, irawdata):
                    # sanity check 2: a pair of chromosomes has to be from a single
                    # individual.
                    if v1[2] != v2[2]:
                        raise ValueError("Chromosomes come from different individual")
                    ids.append(v1[2])
                    inbgens.append(int(float(v1[3])))
                    genos.append([[g1, g2] for g1, g2 in zip(v1[5:], v2[5:])])

                samples.append(FullSample(fname, ids, genos, inbgens))

        return samples

    @staticmethod
    def fromjson(fname):
        """
        Construct an instance of FullSample from a JSON file.  This file should
        have already stored all the necessary information and nothing more.

        This function always return a single FullSample object.
        """

        with open(fname, "r") as f:
            data = json.load(f)
            ids = [val[0] for val in data]
            inbgens = [val[1] for val in data]
            genos = [val[2] for val in data]
            s = FullSample(fname, ids, genos, inbgens)
        return s

    def tojson():
        """
        Write this sample to json-formatted string.
        """

        data = [[i, j, k] for i, j, k in zip(self._ids, self._inbgens, self._genos)]

        return json.dumps(data)

def tonexus(ss, miss, sep):
    """
    Create a nexus-formatted string from a list of samples.
    """

    if type(ss) is not list:
        ss = [ss]

    npops = len(ss)

    lines = [
            "#nexus",
            "begin gdadata",
            "dimensions npops={} nloci={};".format(npops, ss[0].nloc),
            "format missing={} separator={};".format(miss, sep),
            "matrix"
            ]

    for pop, s in enumerate(ss):
        lines.append("{};".format(s.source))
        for idx, geno in zip(s.ids, s.genotypes):
            line = [str(idx)] + ["{}{}{}".format(i[0], sep, i[1]) for i in geno]
            lines.append(" ".join(line))
        if pop < npops - 1:
            lines.append(",")
        else:
            lines.append(";");

    lines.append("end;")

    return lines

def tormes(ss):
    """
    Create a RMES-formatted string from a list of samples.

    The output can be fed to RMES.
    """

    if type(ss) is not list:
        ss = [ss]

    npops = len(ss)

    lines = [str(npops)]

    for i, s in enumerate(ss):
        lines.extend([str(s.source), str(s.n), str(s.nloc)])

    for s in ss:
        for g in s.genotypes:
            lines.append(" ".join([i[0] != i[1] for i in g]))

    return lines

