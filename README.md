# Firedrake-Fluids

A collection of numerical models for fluid flow simulation, using the [Firedrake](http://www.firedrakeproject.org) performance-portable automated solution framework.

[![Documentation Status](https://readthedocs.org/projects/firedrake-fluids/badge/?version=latest)](https://readthedocs.org/projects/firedrake-fluids/?badge=latest)

## Dependencies

Firedrake-Fluids depends on:
* [Firedrake](http://firedrakeproject.org)
* firedrake-adjoint (part of the [dolfin-adjoint](http://dolfin-adjoint.org/) package)
* [libspud](https://launchpad.net/spud)
* [pytest](http://pytest.org) (optional, but required to run `make test`)
* [gmsh](http://geuz.org/gmsh/) (optional, but required to run some tests)
* [Sphinx](http://sphinx-doc.org/) (optional, but required to build the documentation from source)

Note that the codebase is designed to run on the Linux operating system. All development and testing of Firedrake-Fluids is done on the Ubuntu Linux distribution.

## Quick start
* Install the Firedrake-Fluids Python module using

```
sudo python setup.py install
```

Alternatively, append the `firedrake_fluids` directory to your `PYTHONPATH` environment variable using e.g.:

```
export PYTHONPATH=$PYTHONPATH:/home/christian/firedrake-fluids/firedrake_fluids
```

* Run the tests (which are located in the `tests` directory) using

```
make test
```

from the Firedrake-Fluids base directory, to ensure that Firedrake-Fluids is working properly.

* You can setup a new shallow water simulation using Diamond, which comes with the libspud dependency.

```
diamond -s schemas/shallow_water.rng my_new_simulation_setup_file.swml
```

* Run the simulation using

```
python firedrake_fluids/shallow_water.py my_new_simulation_setup_file.swml
```

## Documentation

The Firedrake-Fluids documentation can be compiled using [Sphinx](http://sphinx-doc.org/) by running

```
make html
```

from the Firedrake-Fluids `docs` directory. Alternatively, the latest version of the documentation can be viewed [here](http://firedrake-fluids.readthedocs.org/en/latest/).

## Citing

To cite Firedrake-Fluids, please use:

**C. T. Jacobs and M. D. Piggott (2015)**. *Firedrake-Fluids v0.1: numerical modelling of shallow water flows using an automated solution framework*, Geoscientific Model Development 8(3):533-547, DOI:[10.5194/gmd-8-533-2015](http://dx.doi.org/10.5194/gmd-8-533-2015).

## Build status

[Buildbot](http://buildbot.net/) is used to perform automated testing of Firedrake-Fluids each time a change is made to the codebase. The current status of the `master` branch can be found [here](http://buildbot.ese.ic.ac.uk:8080/builders/trusty-firedrake-fluids/).

## Contact

If you have any questions about Firedrake-Fluids, please send an email to <c.jacobs10@imperial.ac.uk>.

## License

Firedrake-Fluids is released under the GNU General Public License. Please see the file called `COPYING` for more information.

