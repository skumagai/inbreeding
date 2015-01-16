# -*- mode: python; coding: utf-8; -*-

# test_record_selfing_infinite_sites.py - Tests for recording the generations of
# continuous selfing before the first outcrossing event under the
# infinite sites model.

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

import simuOpt
simuOpt.setOptions(quiet=True, alleleType='binary')
import simuPOP as simu

import partial_selfing.common as cf
import partial_selfing.infinite_sites as isf

class TestRecordSelfing:

    def setUp(self):
        # A locus with 10 sites
        "Let there are 5 loci with 2 sites each"
        self.allele_length = 2
        self.loci = 5
        self.pop = simu.Population(size=10,
                                   loci=self.allele_length * self.loci,
                                   infoFields='self_gen')
        self.sim = simu.Simulator(pops = self.pop)
        self.initOps = [simu.InitSex(sex=[simu.MALE, simu.FEMALE]),
                        simu.InitInfo(0, infoFields=['self_gen'])]
        self.sexMode = (simu.GLOBAL_SEQUENCE_OF_SEX, simu.MALE, simu.FEMALE)

    def test_pure_outcrossing(self):
        """Test pure outcrossing population.

        Without any selfing in history, the infoField `self_gen` should have
        0 for all loci for all individuals.
        """
        self.sim.evolve(
            initOps = self.initOps,
            matingScheme=simu.RandomMating(subPopSize = 10,
                                           sexMode = self.sexMode,
                                           ops = [
                                               simu.MendelianGenoTransmitter(),
                                               cf.MyOutcrossingTagger()
                                           ]),
            gen = 10)

        for pop in self.sim.populations():
            for ind in pop.individuals():
                assert ind.info('self_gen') == 0

    def test_pure_selfing(self):
        """Test pure inbreeding population.

        When all individuals undergo selfing, the values in `self_gen` should be identical
        to the number of generations since simulations started.
        """
        self.sim.evolve(
            initOps = self.initOps,
            matingScheme = simu.SelfMating(subPopSize = 10,
                                           sexMode = self.sexMode,
                                           ops = [
                                               simu.SelfingGenoTransmitter(),
                                               cf.MySelfingTagger()
                                               ]),
            gen = 10)

        for pop in self.sim.populations():
            for ind in pop.individuals():
                assert ind.info('self_gen') == 10


    def test_selfing_after_outcrossing(self):
        """Test scinario: generations of selfing after generations of outcrossing."""

        # outcrossing
        self.sim.evolve(
            initOps = self.initOps,
            matingScheme=simu.RandomMating(subPopSize = 10,
                                           sexMode = self.sexMode,
                                           ops = [
                                               simu.MendelianGenoTransmitter(),
                                               cf.MyOutcrossingTagger()
                                           ]),
            gen = 10)

        for pop in self.sim.populations():
            for ind in pop.individuals():
                assert ind.info('self_gen') == 0

        assert pop.dvars().gen == 10

        # selfing
        self.sim.evolve(
            initOps = self.initOps,
            matingScheme = simu.SelfMating(subPopSize = 10,
                                           sexMode = self.sexMode,
                                           ops = [
                                               simu.SelfingGenoTransmitter(),
                                               cf.MySelfingTagger()
                                               ]),
            gen = 10)


        for pop in self.sim.populations():
            for ind in pop.individuals():
                assert ind.info('self_gen') == 10

        assert pop.dvars().gen == 20


    def test_outcrossing_after_selfing(self):
        """Test scinario: generations of outcrossing after generations of selfing."""

        # selfing
        self.sim.evolve(
            initOps = self.initOps,
            matingScheme = simu.SelfMating(subPopSize = 10,
                                           sexMode = self.sexMode,
                                           ops = [
                                               simu.SelfingGenoTransmitter(),
                                               cf.MySelfingTagger()
                                               ]),
            gen = 10)

        for pop in self.sim.populations():
            for ind in pop.individuals():
                assert ind.info('self_gen') == 10

        assert pop.dvars().gen == 10

        # outcrossing
        self.sim.evolve(
            initOps = self.initOps,
            matingScheme=simu.RandomMating(subPopSize = 10,
                                           sexMode = self.sexMode,
                                           ops = [
                                               simu.MendelianGenoTransmitter(),
                                               cf.MyOutcrossingTagger()
                                           ]),
            gen = 10)

        for pop in self.sim.populations():
            for ind in pop.individuals():
                assert ind.info('self_gen') == 0

        assert pop.dvars().gen == 20


    def test_mixed_mating_pure_outcrossing(self):
        """
        Test scenario: mixture mating but selfing rate is set to zero.

        The population is effectively pure outcrossing.
        """
        ms = isf.get_mating_operator(0, 0, 10, loci=self.loci, allele_length=self.allele_length)

        self.sim.evolve(
            initOps = self.initOps,
            matingScheme=ms,
            gen = 10)

        for pop in self.sim.populations():
            for ind in pop.individuals():
                assert ind.info('self_gen') == 0


    def test_mixed_mating_pure_selfing(self):
        """
        Test scenario: mixture mating but selfing rate is set to one.

        The population is effectively pure outcrossing.
        """
        ms = isf.get_mating_operator(0, 1, 10, loci=self.loci, allele_length=self.allele_length)

        self.sim.evolve(
            initOps = self.initOps,
            matingScheme=ms,
            gen = 10)

        for pop in self.sim.populations():
            for ind in pop.individuals():
                assert ind.info('self_gen') == 10


    # TODO: find a way to test for really mixed mating (0 < selfing
    # rate < 1).
