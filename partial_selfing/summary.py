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
# 2. [ ] number of alleles per locus
# 3. [ ] heterozygosity
# 4. [ ] number of alleles per locus given generation of selfing
# 5. [ ] heterozygosity given generation of selfing
# 6. [ ] generation of selfing given number of alleles
# 7. [ ] generation of selfing given heterozygosity


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

field_types = ('a18','i8','i8','i8','i8','f4','f4','f4','i8','f4','i8','i8')

def selfing_gen(data):
    cols = [col for col in data.columns
            if col[:5] == "locus"]
    # print (cols)
    grouped = data.groupby(rep_keys)

    d = np.zeros((len(grouped),), dtype=[(i,j) for i, j in zip(rep_keys, field_types)])
    i = 0
    for name, group in grouped:
        d[i] = name
        i += 1

    base = pd.DataFrame(d, columns=rep_keys, index=range(len(grouped)))

    gens = [group[group["chromosome"] == 0]["number of selfing"].value_counts()
            for name, group in grouped]
    gens = pd.concat(gens, axis=1).fillna(0).T

    print pd.concat([base, gens], axis=1)
    print(grouped.aggregate("count"))

    for name, group in grouped:
        # print(group.loc[:,rep_keys])
        first = group[0::2]
        del first['chromosome']
        second = group[1::2]
        del second['chromosome']
        second = second.reindex_like(first)
        match = (first != second)

    gens = [group["number of selfing"][::2].value_counts()
            for name, group in grouped]

if __name__ == '__main__':
    import os.path

    file_base = os.path.join("kmar_sim", "theta_0.2_s_0.2.{}.csv")
    files = [file_base.format(i) for i in [1] + range(10,20)]
    # data = [pd.read_csv(f) for f in files]
    # d = pd.concat(data)
    for f in files[0:1]:
        data = pd.read_csv(f)
        selfing_gen(data)
