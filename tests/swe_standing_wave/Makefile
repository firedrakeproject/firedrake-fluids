input:	clean
	gmsh -2 src/swe_standing_wave.geo -o swe_standing_wave.msh

run:
	python ../../models/shallow_water.py swe_standing_wave.swml

clean:
	rm -f *.pvd *.vtu *.stat *.log *~ *.pyc nohup.out *.node *.ele *.edge
