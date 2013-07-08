# -*- mode: python; coding: utf-8; -*-

# db.py - store simulation results in SQLite database

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


import sqlite3

def create_db(db):
    """Create a new database and load schema for storing simulation results."""
    conn = sqlite3.connect(db)
    c = conn.cursor()

    # create tables
    c.execute('''
    -- Store infinite-alleles (IA) genotypes
    CREATE TABLE IF NOT EXISTS ia_genotype(
           ind_id INTEGER REFERENCES individual ON DELETE CASCADE,
           locus INTEGER NOT NULL,
           chromosome INTEGER NOT NULL,
           genotype INTEGER NOT NULL,
           UNIQUE(ind_id, locus, chromosome));

    -- Store infinite-sites (IS) genotypes
    CREATE TABLE IF NOT EXISTS is_genotype(
           ind_id INTEGER REFERENCES individual ON DELETE CASCADE,
           locus INTEGER,
           chromosome INTEGER,
           genotype TEXT,
           UNIQUE(ind_id, locus, chromosome));

    -- Store association between individuals and a run of simulations
    CREATE TABLE IF NOT EXISTS individual(
           id INTEGER PRIMARY KEY ASC,
           run_id INTEGER REFERENCES run ON DELETE CASCADE);

    -- Store association of simulation runs and their parameters.
    CREATE TABLE IF NOT EXISTS run(
           id INTEGER PRIMARY KEY ASC,
           param_id INTEGER REFERENCES param ON DELETE CASCADE);

    -- Store parameters of simulation
    CREATE TABLE IF NOT EXISTS param(
           id INTEGER PRIMARY KEY ASC,
           num_loci INTEGER NOT NULL, -- number of loci,
           pop_size INTEGER NOT NULL, -- population size
           mut_rate REAL NOT NULL, -- mutation rate
           rec_rate REAL NOT NULL, -- recombination rate
           mut_model TEXT CHECK (mut_model IN ("IA", "IS")), -- mutational model (IA or IS)
    );
    ''')
    c.close()
