input:	clean
	./copy_script
run:
	python ../../firedrake_fluids/shallow_water.py MMS_A.swml
	#python ../../firedrake_fluids/shallow_water.py MMS_B.swml
	#python ../../firedrake_fluids/shallow_water.py MMS_C.swml

clean:
	rm -f *.pvd *.vtu *.stat *_A* *_B* *_C* *.log matrix* *~ *.pyc nohup.out *.node *.ele *.edge
