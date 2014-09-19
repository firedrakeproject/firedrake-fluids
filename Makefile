#!/bin/sh

.PHONY: docs

input: 	clean docs

docs:
	@echo ***Compiling the documentation
	cd docs; make html; cd ..

test:
	@echo ***Running the tests
	python -m pytest
	
clean:
	@echo ***Cleaning the docs directory
	cd docs; make clean; cd ..

