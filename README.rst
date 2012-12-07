Simulation of Partial-Selfing (Plant) Populations
=================================================

inbreeding.py
-------------

Conduct forward simulation of partial-selfing population at two linked loci.
This is implemented with simuPOP_, and tested with version 1.0.8.

.. _simuPOP: http://simupop.sourceforge.net/

   usage: partialSelfing.py [-h] [-r R] [-s S] [-m [M [M ...]]] [-b B] [-n NLOCI]
                         [--num-segre-sites NSITES] [--seed SEED] [--explore]
                         POP NGEN NREP [OUTPUT]

   run partial selfing simulations

   positional arguments:
     POP                   population size
     NGEN                  number of generations excluding burn-in period
     NREP                  number of replicates
     OUTPUT                output file name (default: STDOUT)

   optional arguments:
     -h, --help            show this help message and exit
     -r R, --recombination-rate R
                           recombination rate (default: 0.5)
     -s S, --selfing-rate S
                           selfing rate (default: 0)
     -m [M [M ...]], --mutation-rate [M [M ...]]
                           mutation rate (default:0)
     -b B, --burn-in B     burn-in (default: 0)
     -n NLOCI, --num_loci NLOCI
                           number of loci (default: 2)
     --num-segre-sites NSITES
                           maximum number of segregating sites per locus
                           (default: 256)
     --seed SEED           random number seed (default: use posix time)
     --explore             record heterozygosity and number of segregating sites
                           each generation for later inspection to determine an
                           appropriate durtion of burn-in
     --infinite-alleles    use the infinite-alleles model instead of the
                           infinite-sites model

Population heterozygosities for each run are printed to STDOUTAll in CSV format.

Assumptions of this simulation is:

* There is no population structure.
* There is no natural-selection.
* The population is hermaphroditic.
* Population size is constant at ``pop_size``.
* Mutation rate at the first locus is ``mut_rate0``, and at the second locus is ``mut_rate1``.
* A fraction (``selfing_rate``) of organisms are selfers, and they fertilize by themselves.
* The rest (``1 - self_rate``) are outcrossers, and they cannot self-fertilize.
* Identity of parent(s) does not determine if an organism is a selfer or outcrosser.
* Recombination rate between two loci is ``recomb_rate``.
* Evolution is simulated for ``ngen`` generations.
* Simulations are repeated for ``nrep`` times.
