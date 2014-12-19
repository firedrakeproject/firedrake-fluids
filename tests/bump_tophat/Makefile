input:	clean preprocess

preprocess:
	gmsh -2 src/square.geo -o src/square.msh

run:
	python ../../firedrake_fluids/shallow_water.py bump.swml
	python ../../firedrake_fluids/shallow_water.py tophat.swml

clean:
	rm -f *.pvd *.vtu *.stat *_A* *_B* *_C* *.log matrix* *~ *.pyc nohup.out *.node *.ele *.edge
