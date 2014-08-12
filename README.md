# Firedrake-Fluids

A collection of numerical models for fluid flow simulation using the [Firedrake](http://www.firedrakeproject.org) framework.

## Dependencies

Firedrake-Fluids depends on:
* [Firedrake](http://firedrakeproject.org)
* [libspud](https://launchpad.net/spud)
* [pytest](http://pytest.org) (optional, but required to run `make tests`)
* [gmsh](http://geuz.org/gmsh/) (optional, but required to run some tests)
* pdfLaTeX (optional, but required to build the user manual from source)

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

* You can setup a new shallow water simulation using Diamond, which comes with the libspud install.

```
diamond -s schemas/shallow_water.rng my_new_simulation_setup_file.swml
```

* Run the simulation using

```
python models/shallow_water.py my_new_simulation_setup_file.swml
```

## Documentation

The Firedrake-Fluids user manual is written in LaTeX. The .tex source file can be compiled by running

```
make manual
```

from the Firedrake-Fluids base directory.

## Contact

If you have any questions about Firedrake-Fluids, please send an email to <c.jacobs10@imperial.ac.uk>.

## License

Firedrake-Fluids is released under the GNU General Public License. Please see the file called COPYING for more information.

