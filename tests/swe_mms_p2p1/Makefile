input:	clean
	./copy_script
run:
	python ../../models/shallow_water.py MMS_A.swml
	#python ../../models/shallow_water.py MMS_B.swml
	#python ../../models/shallow_water.py MMS_C.swml

clean:
	rm -f *.pvd *.vtu *.stat *_A* *_B* *_C* *_D* src/*A.* src/*B.* src/*C.* src/*D.* *.log matrix* *~ *.pyc nohup.out *.node *.ele *.edge
