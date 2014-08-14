#!/bin/sh

input: 	clean manual

manual:
	@echo **********Compiling the user manual
	cd doc; pdflatex manual.tex; bibtex manual; pdflatex manual.tex; pdflatex manual.tex; cd ..

test:
	@echo **********Running the tests
	python -m pytest
	
clean:
	@echo **********Cleaning doc directory
	cd doc; rm -rf *.log *.aux *.dvi *.pdf *.ps *.toc *.out *.bbl *.fdb_latexmk; cd ..

