input:	clean
	./copy_script
	gmsh -2 src/MMS_A.geo
	gmsh -2 src/MMS_B.geo
	gmsh -2 src/MMS_C.geo

run:
	python ../../firedrake_fluids/shallow_water.py MMS_A.swml
	python ../../firedrake_fluids/shallow_water.py MMS_B.swml
	python ../../firedrake_fluids/shallow_water.py MMS_C.swml

clean:
	rm -f *.vtu *.stat *.s *_A* *_B* *_C* *_D* src/*A.* src/*B.* src/*C.* src/*D.* *.log *~ *.pyc nohup.out *.node *.ele *.edge
