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


import glob
import pandas as pd
import numpy as np

index = ["mutation model",
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
         "generation",
         "individual",
         "chromosome"]

summary_index = ["mutation model",
                 "number of individuals",
                 "number of generations",
                 "number of replicates",
                 "number of loci",
                 "mutation rate",
                 "selfing rate",
                 "recombination rate",
                 "number of burnin generations",
                 "generation"]

file_names = ["selfing",
              "alleles",
              "h_haplo",
              "a_given_s",
              "h_given_s",
              "s_given_h"]

def run(args):
    files = glob.glob(args.PATTERN)
    dfs = process_single_file(files[0])
    for f in files[1:]:
        dfs_tmp = process_single_file(f)
        dfs = [pd.concat([df, df_tmp]) for df, df_tmp in zip(dfs, dfs_tmp)]

    dfs = [df.fillna(0) for df in dfs]
    for df, f in zip(dfs, file_names):
        df.to_csv(f + ".csv", cols=df.columns)

    df = dfs[0].groupby(level=summary_index + ['number of selfing']).sum()
    df.to_csv(file_names[0] + ".summary.csv", cols = df.columns)
    df = dfs[1].groupby(level=summary_index).agg([np.mean, np.std])
    df.to_csv(file_names[1] + ".summary.csv", cols = df.columns)
    df = dfs[2].groupby(level=summary_index + ['heterozygosity']).sum()
    df.to_csv(file_names[2] + ".summary.csv", cols = df.columns)
    df = dfs[3].groupby(level=summary_index + ['number of selfing']).agg([np.mean, np.std])
    df.to_csv(file_names[3] + ".summary.csv", cols = df.columns)
    df = dfs[4].groupby(level=summary_index + ['number of selfing', 'heterozygosity']).sum()
    df.to_csv(file_names[4] + ".summary.csv", cols = df.columns)
    df = dfs[5].groupby(level=summary_index + ['heterozygosity', 'number of selfing']).sum()
    df.to_csv(file_names[5] + ".summary.csv", cols = df.columns)


def process_single_file(f):
    data = pd.read_csv(f, index_col = index)
    return (selfing_gen(data),
            number_of_alleles(data),
            heterozygosity(data),
            number_of_alleles_given_selfing(data),
            heterozygosity_given_selfing(data),
            selfing_given_heterozygosity(data))



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

    return pd.DataFrame(d, columns=keys, index=range(len(grouped))).fillna(0)


def selfing_gen(data):
    d = data[0::2].groupby(level=index[:-2])['number of selfing'].value_counts()

    d.index.names = d.index.names[:-1] + ['number of selfing']
    return pd.DataFrame({'count': d})


def number_of_alleles(data):
    loci = [col for col in data.columns if col[:5] == "locus"]

    return data.groupby(level=index[:-2])[loci].agg(pd.Series.nunique).fillna(0)


def heterozygosity(data):
    loci = [col for col in data.columns if col[:5] == "locus"]
    d = (data[0::2].reset_index(level=['chromosome'])[loci] !=
            data[1::2].reset_index(level=['chromosome'])[loci]) \
        .fillna(0) \
        .sum(axis=1) \
        .groupby(level=index[:-2]) \
        .value_counts()

    d.index.names = d.index.names[:-1] + ['heterozygosity']
    return pd.DataFrame({'count': d})


def number_of_alleles_given_selfing(data):
    loci = [col for col in data.columns if col[:5] == "locus"]

    d = data.set_index('number of selfing', append=True) \
            .reset_index(level=['individual', 'chromosome']) \
            .sortlevel()
    return d.groupby(level=index[:-2] + ['number of selfing'])[loci] \
            .agg(pd.Series.nunique) \
            .fillna(0)


def heterozygosity_given_selfing(data):
    loci = [col for col in data.columns if col[:5] == "locus"]

    d = data.set_index('number of selfing', append=True) \
            .reset_index(level=['individual', 'chromosome']) \
            .sortlevel()[loci]

    d = (d[0::2] != d[1::2]) \
        .fillna(0) \
        .sum(axis=1) \
        .groupby(level=index[:-2] + ['number of selfing']) \
        .value_counts()

    d.index.names = d.index.names[:-1] + ['heterozygosity']
    return pd.DataFrame({'count': d})


def selfing_given_heterozygosity(data):
    # Create a new data frame with information about heterozygosity
    # and without genotype
    loci = [col for col in data.columns if col[:5] == "locus"]

    d = data.set_index('number of selfing', append=True) \
            .reset_index(level=['individual', 'chromosome']) \
            .sortlevel()[loci]

    d = (d[0::2] != d[1::2]) \
        .fillna(0) \
        .sum(axis=1) \
        .reset_index(level=['number of selfing'])
    d.columns = ['number of selfing', 'heterozygosity']

    d = d.set_index('heterozygosity', append=True) \
         .sortlevel()['number of selfing'] \
         .groupby(level=index[:-2] + ['heterozygosity']) \
         .value_counts()

    d.index.names = d.index.names[:-1] + ['number of selfing']
    return pd.DataFrame({'count': d})


if __name__ == '__main__':
    import os.path

    file_base = os.path.join("kmar_sim", "theta_0.2_s_0.2.{}.csv")
    files = [file_base.format(i) for i in [1] + range(10,20)]
    for f in files[0:1]:
    # for f in files:
        data = pd.read_csv(f, index_col = index)
        print selfing_gen(data).head()
        # print number_of_alleles(data).head()
        # print heterozygosity(data).head()
        # print number_of_alleles_given_selfing(data).head()
        # print heterozygosity_given_selfing(data).head()
        # print selfing_given_heterozygosity(data).head()
