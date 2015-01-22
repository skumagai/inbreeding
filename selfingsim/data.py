"""
selfingsim.data
===============

Representation of samples and its associated methods.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

try:
    str = unicode
except NameError:
    pass

# standeard imports
from collections import Counter
import csv
import io
from itertools import groupby, izip
import json
import os.path
import random

def createsample(fname, gen=None):
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
        """
        Returns the name of a file containing this sample.
        """
        return self._src

    @property
    def ids(self):
        """
        Returns IDs in the sample.
        """
        return self._ids

    @property
    def genotypes(self):
        """
        Returns a list of genotypes in the sample.
        """
        return self._genos

    @property
    def nsam(self):
        """
        Returns the sample size.
        """
        return self._nsam

    @property
    def nloc(self):
        """
        Returns the number of loci in the sample.
        """
        return self._nloc

    def sample(self, nsam, locs=None):
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
            nloc = len(self._genos[0])
            if len([True for i in locs if i >= nloc]) > 0:
                raise ValueError("Specified locus outside of range")
            ret = locs
        return ret


    def _drawindividuals(self, newnsam):
        """
        Draw n individuals without replacement from the current sample.
        """
        return random.sample(xrange(self._nsam), newnsam)

    @staticmethod
    def fromphase(fname):
        """
        Create an instance of BasicSample.
        """
        with io.open(fname, "rU") as fhandle:
            next(fhandle)
            nloc = int(next(fhandle).strip())
            types = next(fhandle).strip()

            # Sanity check
            if len(types) != nloc:
                raise ValueError("Inconsistent number of loci.")

            ids = []
            genos = []

            for line in fhandle:
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

    def tonexus(self, miss="?", sep="/"):
        """
        Write this sample as a nexus-formatted string.
        """
        return tonexus([self], miss, sep)

    def tormes(self):
        """
        Write this sample as a RMES-formatted string.
        """
        return tormes([self])

    def _afreqs(self):
        """
        Return allele frequency spectrum per-locus.
        """
        alleles = [[gene for ind in self._genos for gene in ind[loc]] for loc in xrange(self._nloc)]
        return [[v / (2 * self._nsam) for v in Counter(a).values()] for a in alleles]

    def _hobs(self):
        """
        Obtains per-locus observed heterozygosity.
        """
        genos = [[ind[loc] for ind in self._genos] for loc in xrange(self._nloc)]
        return [sum(1 if g[0] != g[1] else 0 for g in geno) / self.nsam for geno in genos]

    @staticmethod
    def _hexp(afreqs):
        """
        Computes per-locus expected heterozygosity from allele frequency.
        """
        return [1 - sum(a * a for a in locus) for locus in afreqs]

    def _fis(self, hobs, hexp, correct=False):
        """
        Computes inbreeding coefficient (Fis) per locus.
        """
        ngenes = 2 * self._nsam
        if correct:
            fis = [(exp - obs + obs / ngenes) / (exp - obs / ngenes)
                   if exp - obs / ngenes != 0. else float('nan')
                   for obs, exp in zip(hobs, hexp)]
        else:
            fis = [1 - obs / exp if exp != 0. else float('nan') for obs, exp in zip(hobs, hexp)]

        return fis

    def inbreedingcoefficient(self):
        """
        Compute per-locus and overall inbreeding coefficient Fis.

        In addition to Fis, this method also returns other related statistics
        (observed and expected heterozygosities, and bias-corrected Fis as well as
        number of observed alleles.
        """
        afreqs = self._afreqs()
        hobs = self._hobs()
        hexp = self._hexp(afreqs)

        nalleles = [len(a) for a in afreqs]

        fis = self._fis(hobs, hexp, False)
        fis2 = self._fis(hobs, hexp, True)

        fis.append(self._fis([sum(hobs)], [sum(hexp)], False)[0])
        fis2.append(self._fis([sum(hobs)], [sum(hexp)], True)[0])

        hobs.append(sum(hobs) / len(hobs))
        hexp.append(sum(hexp) / len(hexp))
        nalleles.append(sum(nalleles) / len(nalleles))

        keys = [i for i in xrange(self.nloc)]
        keys.append("overall")

        data = []
        for vals in zip(keys, hobs, hexp, fis, fis2, nalleles):
            data.append({
                "key": vals[0],
                "hobs": vals[1],
                "hexp": vals[2],
                "fis": vals[3],
                "fisc": vals[4],
                "nalleles": vals[5]})

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
        """
        Returns generation since the last outcrossing per-individual.
        """
        return self._inbgens

    def sample(self, nsam, locs=None):
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
        with io.open(fname, "r") as fhandle:
            reader = csv.reader(fhandle, delimiter=str("\t"))
            # Throw out a header row
            next(fhandle)
            rows = [row for row in reader if int(row[1]) == gen]

            # sanity check 1: total number of chromosomes (rows) must be multiple
            # of 2, as diploids have two chromosomes.
            if len(rows) % 2 != 0:
                raise ValueError("Number of chromosomes not mulitple of 2")


            for _, rawdata in groupby(rows, lambda x: x[0]):
                irawdata = iter(rawdata)

                ids = []
                inbgens = []
                genos = []
                for vals in izip(irawdata, irawdata):
                    # sanity check 2: a pair of chromosomes has to be from a single
                    # individual.
                    if vals[0][2] != vals[1][2]:
                        raise ValueError("Chromosomes come from different individual")
                    ids.append(vals[0][2])
                    inbgens.append(int(float(vals[0][3])))
                    genos.append([[gvals[0], gvals[1]] for gvals in zip(vals[0][5:], vals[1][5:])])

                samples.append(FullSample(fname, ids, genos, inbgens))

        return samples

    @staticmethod
    def fromjson(fname):
        """
        Construct an instance of FullSample from a JSON file.  This file should
        have already stored all the necessary information and nothing more.

        This function always return a single FullSample object.
        """

        with io.open(fname, "r") as fhandle:
            data = json.load(fhandle)
        ids = [val[0] for val in data]
        inbgens = [val[1] for val in data]
        genos = [val[2] for val in data]
        return FullSample(fname, ids, genos, inbgens)

    def tojson(self):
        """
        Write this sample to json-formatted string.
        """

        data = [[i, j, k] for i, j, k in zip(self._ids, self._inbgens, self._genos)]

        return str(json.dumps(data))

def tonexus(samples, miss, sep):
    """
    Create a nexus-formatted string from a list of samples.
    """

    if type(samples) is not list:
        samples = [samples]

    npops = len(samples)

    lines = [
        "#nexus",
        "begin gdadata",
        "dimensions npops={} nloci={};".format(npops, samples[0].nloc),
        "format missing={} separator={};".format(miss, sep),
        "matrix"]

    for pop, sample in enumerate(samples):
        lines.append("{}:".format(sample.source))
        for idx, geno in zip(sample.ids, sample.genotypes):
            line = [str(idx)] + ["{}{}{}".format(i[0], sep, i[1]) for i in geno]
            lines.append(" ".join(line))
        if pop < npops - 1:
            lines.append(",")
        else:
            lines.append(";")

    lines.append("end;")

    return lines

def tormes(samples):
    """
    Create a RMES-formatted string from a list of samples.

    The output can be fed to RMES.
    """

    if type(samples) is not list:
        samples = [samples]

    npops = len(samples)

    lines = [str(npops)]

    for sample in samples:
        lines.extend([str(sample.source), str(sample.nsam), str(sample.nloc)])

    for sample in samples:
        for genotype in sample.genotypes:
            lines.append(" ".join(["1" if geno[0] != geno[1] else "0" for geno in genotype]))

    return lines

