"""
selfingsim.simulate
###################

Run a forward-in-time population genetics simulation for partially selfing organisms.
"""


# standard imports
import argparse
import json
import sys


def run():
    """
    Runs simulations as a stand-alone script.
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    setup_command_line(subparsers)
    args = parser.parse_args()
    args.func(args)


def setup_command_line(subparsers):
    """
    Sets up command line interface.
    """
    # Forward simulation
    parser = subparsers.add_parser(
        "simulate", help="run forward simulations with selfing"
    )
    parser.add_argument(
        "config",
        type=argparse.FileType("r"),
        help="the settings of simulation as a JSON file",
    )
    parser.add_argument(
        "substs",
        type=str,
        nargs="*",
        help="substitutions plugged into an output file name (specified in config)",
    )
    parser.set_defaults(func=simulate)


def simulate(args):
    """
    Runs simulation.
    """
    config = Config(json.load(args.config), args.substs)

    if config.mutation_model == "infinite sites":
        print("The infinite-sites model is disabled")
        sys.exit(1)
        exec_infinite_sites(config)
    elif config.mutation_model == "infinite alleles":
        exec_infinite_alleles(config)
    else:
        print("Unknown mutational model specified", file=sys.stderr)
        sys.exit(1)


class Config(object):
    """
    Stores settings of simulations.
    """

    def __init__(self, cobj, substs):
        self._params = {}
        # Sets up simple parameters
        self._addparam(cobj, "population", "N")
        self._addparam(cobj, "population", "loci")
        self._addparam(cobj, "population", "r")
        self._addparam(cobj, "general", "outfile", lambda x: x.format(*substs))
        self._addparam(cobj, "general", "gens", lambda x: self._params["N"] * x)
        self._addparam(cobj, "general", "burnin", lambda x: self._params["N"] * x)
        self._addparam(cobj, "general", "debug")

        # check if "output per" exists in an input file.  If not, set the value to
        # the last generation.
        self._params["output_per"] = self._params["N"]
        try:
            self._params["output_per"] *= cobj["general"]["output per"]
        except KeyError:
            self._params["output_per"] *= self._params["gens"] * self._params["burnin"]

        # Sets up more complex parameters
        # start with mtaing scheme
        self._addmating(cobj)
        # then mutation model
        self._addmutation(cobj)
        # finally, initial genotpye
        self._addinitgenotype(cobj)

    def _addinitgenotype(self, cobj):
        """
        Adds settings of initial genotypes of a simulated population.
        """
        init = cobj["population"]["init"]
        if init == "unique" or init == "monomorphic":
            self._params["initial_genotype"] = [init]
        elif type(init) is int and init > 0:
            self._params["initial_genotype"] = ["count", init]
        else:
            try:
                if all(i >= 0.0 for i in init):
                    norm = sum(init)
                    self._params["initial_genotype"] = [
                        "frequency",
                        [i / norm for i in init],
                    ]
                else:
                    sys.exit(
                        "Initial frequencies of alleles must be all non-negative: {}.".format(
                            init
                        )
                    )
            except TypeError:
                sys.exit("Unknown init: {}.".format(init))

    def _addparam(self, cobj, sec, key, mod=None):
        """
        Adds settings of simple parameters.
        """
        mkey = key.replace(" ", "_")
        try:
            if mod != None:
                self._params[mkey] = mod(cobj[sec][key])
            else:
                self._params[mkey] = cobj[sec][key]
        except KeyError:
            sys.exit('"{}" not found in a config file.'.format(key))

    def _addmutation(self, cobj):
        """
        Adds settings related to mutations.
        """
        try:
            mutation = cobj["population"]["mutation"]
        except KeyError:
            sys.exit("Settings for mutational model not found in a config file.")

        try:
            model = mutation["model"]
            npop = self._params["N"]
            if model == "infinite alleles":
                self._params["mutation_model"] = model
                self._getmutationrate(mutation["theta"], npop)
                self._params["allele_length"] = 1
            elif model == "infinite sites":
                self._params["mutation_model"] = model
                self._getmutationrate(mutation["theta"], npop)
                self._params["allele_length"] = mutation["allele length"]
            else:
                sys.exit("Unrecognized mutation model.")
        except KeyError:
            sys.exit("Unrecognized mutation model.")

    def _getmutationrate(self, rate, npop):
        """
        Un-scales and sets mutation rates.
        """
        if type(rate) is float:
            self._params["m"] = [rate / (4 * npop)] * self._params["loci"]
        elif type(rate) is list and type(rate[0]) is dict:
            self._params["m"] = [
                r["value"] / (4 * npop) for r in rate for _ in range(r["times"])
            ]
        elif type(rate) is list and type(rate[0]) is float:
            if len(rate) == self._params["loci"]:
                self._params["m"] = [r / (4 * npop) for r in rate]
            else:
                sys.exit("Mutation parametrs not fully specified.")
        else:
            sys.exit("Mutation parameters in wrong format.")

    def _addmating(self, cobj):
        """
        Sets mating scheme.
        """
        try:
            mating = cobj["population"]["mating"]
        except KeyError:
            sys.exit("mating scheme not specified.")

        try:
            model = mating["model"]
            self._params["mating_model"] = model
        except KeyError:
            sys.exit("Mating model(scheme) not specified.")

        if model == "pure hermaphroditism":
            self._pure_hermaphroditism(mating)
        elif model == "androdioecy":
            self._androdioecy(mating)
        elif model == "gynodioecy":
            self._gynodioecy(mating)
        else:
            sys.exit('Unrecognized mating model "{}".'.format(model))

    def _pure_hermaphroditism(self, mating):
        try:
            self._params["sstar"] = mating["s*"]
        except KeyError:
            self._params["stilde"] = mating["s tilde"]
            self._params["tau"] = mating["tau"]

    def _androdioecy(self, mating):
        self._pure_hermaphroditism(mating)
        self._params["N_hermaphrodites"] = mating["N_hermaphrodites"]

    def _gynodioecy(self, mating):
        try:
            self._params["sstar"] = mating["s*"]
            self._params["H"] = mating["H"]
        except KeyError:
            self._params["tau"] = mating["tau"]
            self._params["a"] = mating["a"]
            self._params["sigma"] = mating["sigma"]
        self._params["N_hermaphrodites"] = mating["N_hermaphrodites"]

    def __getattr__(self, name):
        """
        Makes self._params[X] accessible as self.X.
        """
        return self._params[name]


# Selectively import simuPOP with an appropriate alleleType.  Because
# python caches imported modules, second or later call of import
# does not simuPOP re-imported.  Therefore, those
# submodules will automatically use the right version of simuPOP.
def exec_infinite_sites(config):
    """
    Launches simulations with the infinite-sites model.
    """
    import selfingsim.infinite_sites as model

    model.run(config)


# See the comment in front of exec_two_loci
def exec_infinite_alleles(config):
    """
    Launches simulations with the infinite-alleles model.
    """
    import selfingsim.infinite_alleles as model

    model.run(config)


if __name__ == "__main__":
    run()
