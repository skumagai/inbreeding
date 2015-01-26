==========
Selfingsim
==========

Selfingsim is a forward-in-time population genetic simulator for studying
evolution of partially selfing organisms.
It was developed for and used in
"Bayesian co-estimation of selfing rate and locus-specific mutation
rates for a partially selfing population"
(Redelings *et al.* 2015, `link`_).

Characteristics of selfingsim are:
- a constant population size
- inter-locus (but not intra-locus) recombination
- three different mating schemes:
  - pure hermaphroditism,
  - androdioecy: hermaphrodites and males,
  - gynodioecy: hermaphrodites and females.
- two mutational models:
  - the infinite-alleles model,
  - the infinite-sites model.

Installation
============

Selfingsim depends on three projects:
1. Python (required; version 2.7 or version 3.0 or later)
2. `simuPOP`_ (rquired): a framework for forward-in-time population genetic simulations
3. `nose`_ (optional): a testing framework for python



## Model of mutations

Selfingsim can simulate a population under two models of mutations:
- Infinite alleles model
- Infinite sites model

Our project used the infinite alleles model, so the code for the model
is well tested.
Each locus under this model is stored in *long* (64-bits), and
a simulation can have at most 2^64-1 mutations from start to end.

The infinite sites model was not used in the project,
and the code for it is only tested very lightly.
Each locus under this model is stored as an array of bits.
Although the length of array is configurable, but dynamic
resizing of the array is not allowed.
A n-length array can store upto n standing mutations.
If all elements are occupied by standing mutations,
the next mutation will make the entire simulator abort.

## Installation

### Dependencies

The simulator is implemented in Python using [simuPOP][],
"a general-purpose individual-based forward-time population genetics
simulation environment".
Unit testing also depends on [nose][].

In order to run selfingsim, simuPOP is always required, whereas nose is
optional.

To install simuPOP, please refer [its download and installation page][simuPOPinst].

## Installation

Run `pip install git+git://github.com/skumagai/skumagai.git` from command line
to install selfingsim.
This method needs you to have installed git.

Alternatively you can download selfingsim from [this link][selfingsimdl], and then
run `pip install master.zip`.

## Running selfingsim

Selfingsim is designed to work from command line.
Basic form of invokation is:

```
selfingsim [command] [options]
```

You can get a list of available commands by `selfingsim -h`.

### Simulation

To run simulations,

```
selfingsim simulate configfile
```

where configfile should be replaced by the name of configuration file
in JSON format.
An example configuration file can be found at [link][selfingsimex].
This example is annotated by ``// comment``, but an actual configuration file
can't have annotation.

.. _link: http://www.example.com
[selfingsimdl]: https://github.com/skumagai/selfingsim/archive/master.zip
[selfingsimex]: https://github.com/skumagai/selfingsim/blob/master/example.json.annotated
.. _simuPOP: http://simupop.sourceforge.net
[simuPOPinst]: http://simupop.sourceforge.net/Main/Download
.. _nose: https://github.com/nose-devs/nose
