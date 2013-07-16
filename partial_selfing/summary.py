# -*- mode: python; coding: utf-8; -*-

# summary.py - summarize simulation result

# Copyright (C) 2013 Seiji Kumagai

# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice (including the next
# paragraph) shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.


# This script summarizes simulation results in per-replicate basis.
# Summary includes:
#
# 1. [X]  generations of selfing
# 2. [X] number of alleles per locus
# 3. [X] heterozygosity
# 4. [X] number of alleles per locus given generation of selfing
# 5. [X] heterozygosity given generation of selfing
# 6. [D] generation of selfing given number of alleles
# 7. [X] generation of selfing given heterozygosity


import pandas as pd
import numpy as np

rep_keys = ["mutation model",
            "number of individuals",
            "number of generations",
            "number of replicates",
            "number of loci",
            "mutation rate",
            "selfing rate",
            "recombination rate",
            "number of burnin generations",
            "random number seed",
            "replicate",
            "generation"]

field_types = ['a18','i8','i8','i8','i8','f4','f4','f4','i8','f4','i8','i8']

def create_base(grouped, add_field=None):
    if add_field is None:
        keys = rep_keys
        fields = field_types
    else:
        keys = rep_keys + [add_field[0]]
        fields = field_types + [add_field[1]]

    d = np.zeros((len(grouped),), dtype=[(i,j) for i, j in zip(keys, fields)])
    i = 0
    for name, group in grouped:
        d[i] = name
        i += 1

    return pd.DataFrame(d, columns=keys, index=range(len(grouped)))

def selfing_gen(data):
    grouped = data.groupby(rep_keys)
    base = create_base(grouped)

    gens = [group[group["chromosome"] == 0]["number of selfing"].value_counts()
            for name, group in grouped]
    gens = pd.concat(gens, axis=1).fillna(0).T

    return pd.concat([base, gens], axis=1)


def number_of_alleles(data):
    loci = [col for col in data.columns if col[:5] == "locus"]

    grouped = data.groupby(rep_keys)
    base = create_base(grouped)

    counts = [group[loci].apply(pd.Series.nunique, axis=0)
              for name, group in grouped]
    counts = pd.concat(counts, axis=1).T

    return pd.concat([base, counts], axis=1)


def heterozygosity(data):
    loci = [col for col in data.columns if col[:5] == "locus"]

    grouped = data.groupby(rep_keys)
    base = create_base(grouped)

    matches = []
    for name, group in grouped:
        first = group[0::2][loci].reset_index(drop=True)
        second = group[1::2][loci].reset_index(drop=True)
        match = (first != second) \
            .apply(pd.Series.value_counts, axis=1) \
            .fillna(0)
        try:
            matches.append(match[True].value_counts())
        except KeyError:
            matches.append(pd.Series([len(first)], index=[0]))

    matches = pd.concat(matches, axis=1).fillna(0).T
    return pd.concat([base, matches], axis=1)


def number_of_alleles_given_selfing(data):
    loci = [col for col in data.columns if col[:5] == "locus"]

    grouped = data.groupby(rep_keys + ["number of selfing"])
    base = create_base(grouped, ["number of selfing", "i8"])

    counts = [group[loci].apply(pd.Series.nunique, axis=0)
              for name, group in grouped]
    counts = pd.concat(counts, axis=1).T.fillna(0)
    return pd.concat([base, counts], axis=1)


def heterozygosity_given_selfing(data):
    loci = [col for col in data.columns if col[:5] == "locus"]

    grouped = data.groupby(rep_keys + ["number of selfing"])
    base = create_base(grouped, ["number of selfing", "i8"])

    matches = []
    for name, group in grouped:
        first = group[0::2][loci].reset_index(drop=True)
        second = group[1::2][loci].reset_index(drop=True)
        match = (first != second) \
            .apply(pd.Series.value_counts, axis=1) \
            .fillna(0)
        try:
            matches.append(match[True].value_counts())
        except KeyError:
            matches.append(pd.Series([len(first)], index=[0]))

    matches = pd.concat(matches, axis=1).fillna(0).T
    return pd.concat([base, matches], axis=1)


def selfing_given_heterozygosity(data):
    # Create a new data frame with information about heterozygosity
    # and without genotype
    loci = [col for col in data.columns if col[:5] == "locus"]

    subdata = data[0::2]
    first = subdata[loci].reset_index(drop=True)
    second = data[1::2][loci].reset_index(drop=True)
    match = (first != second) \
            .apply(pd.Series.value_counts, axis=1) \
            .fillna(0)

    try:
        match = match[True]
    except KeyError:
        match = pd.Series([0] * len(first))

    base = subdata[rep_keys + ['number of selfing']].reset_index(drop=True)
    base['heterozygosity'] = match


    # group by heterozygosity for each replicate separately.
    grouped = base.groupby(rep_keys + ["heterozygosity"])
    b = create_base(grouped, ["heterozygosity", "i8"])

    gens = [group["number of selfing"].value_counts()
            for name, group in grouped]
    gens = pd.concat(gens, axis=1).fillna(0).T

    return pd.concat([b, gens], axis=1)


if __name__ == '__main__':
    import os.path

    file_base = os.path.join("kmar_sim", "theta_0.2_s_0.2.{}.csv")
    files = [file_base.format(i) for i in [1] + range(10,20)]
    for f in files[0:1]:
        data = pd.read_csv(f)
        selfing_gen(data)
        number_of_alleles(data)
        heterozygosity(data)
        number_of_alleles_given_selfing(data)
        heterozygosity_given_selfing(data)
        selfing_given_heterozygosity(data)
