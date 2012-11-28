Simulation of Partial-Selfing (Plant) Populations
=================================================

inbreeding.py
-------------

Conduct forward simulation of partial-selfing population at two linked loci.
This is implemented with simuPOP_, and tested with version 1.0.8.

.. _simuPOP: http://simupop.sourceforge.net/

   usage: inbreeding.sh inbreeding.py pop_size mut_rate0 mut_rate1 selfing_rate recomb_rate ngen nrep

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
