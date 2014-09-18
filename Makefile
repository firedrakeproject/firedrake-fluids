#!/bin/sh

.PHONY: doc

input: 	clean doc

doc:
	@echo ***Compiling the user documentation
	cd doc; make html; cd ..

test:
	@echo ***Running the tests
	python -m pytest
	
clean:
	@echo ***Cleaning doc directory
	cd doc; make clean; cd ..

