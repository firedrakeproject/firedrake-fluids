# Firedrake-Fluids

A collection of numerical models for fluid flow simulation, using the [Firedrake](http://www.firedrakeproject.org) performance-portable automated solution framework.

## Dependencies

Firedrake-Fluids depends on:
* [Firedrake](http://firedrakeproject.org)
* [libspud](https://launchpad.net/spud)
* [pytest](http://pytest.org) (optional, but required to run `make test`)
* [gmsh](http://geuz.org/gmsh/) (optional, but required to run some tests)
* [Sphinx](http://sphinx-doc.org/) (optional, but required to build the user manual from source)

Note that the codebase is designed to run on the Linux operating system. All development and testing of Firedrake-Fluids is done on the Ubuntu Linux distribution.

## Quick start
* Append the Firedrake-Fluids `models` directory to your `PYTHONPATH` environment variable using, for example,

```
export PYTHONPATH=$PYTHONPATH:/home/christian/firedrake-fluids/models
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
python models/shallow_water.py my_new_simulation_setup_file.swml
```

## Documentation

The Firedrake-Fluids documentation can be compiled using [Sphinx](http://sphinx-doc.org/) by running

```
make html
```

from the Firedrake-Fluids `doc` directory.

## Citing
Please cite the Firedrake-Fluids model description paper:

**C. T. Jacobs and M. D. Piggott (2014)**. *Firedrake-Fluids v0.1: numerical modelling of shallow water flows using a performance-portable automated solution framework*, Geosci. Model Dev. Discuss., 7, 5699-5738, [doi:10.5194/gmdd-7-5699-2014](http://dx.doi.org/10.5194/gmdd-7-5699-2014).

## Contact

If you have any questions about Firedrake-Fluids, please send an email to <c.jacobs10@imperial.ac.uk>.

## License

Firedrake-Fluids is released under the GNU General Public License. Please see the file called COPYING for more information.

