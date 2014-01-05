input:	clean
	gmsh -2 -bin src/swe_standing_wave.geo -o swe_standing_wave.msh

run:
	python ../../models/shallow_water.py sw_standing_wave.swml

clean:
	rm -f *.vtu *.stat *.s *.log matrix* *~ *.pyc nohup.out *.node *.ele *.edge
