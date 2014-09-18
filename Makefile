#!/bin/sh

input: 	clean build doc

build:
	@echo "*** Building Firedrake-Fluids"
	sudo python setup.py build > build.log 2>&1 || cat build.log

install:
	@echo "*** Installing Firedrake-Fluids"
	sudo python setup.py install > install.log 2>&1 || cat install.log

doc:
	@echo "*** Compiling the user documentation"
	cd doc; make html; cd ..

test:
	@echo "*** Running the tests"
	python -m pytest
	
clean:
	@echo "*** Cleaning doc directory"
	cd doc; make clean; cd ..

