input:	clean
	./copy_script
run:
	python ../../models/shallow_water.py MMS_A.swml
	#python ../../models/shallow_water.py MMS_B.swml
	#python ../../models/shallow_water.py MMS_C.swml

clean:
	rm -f *.vtu *.stat *.s *_A* *_B* *_C* *_D* src/*A.* src/*B.* src/*C.* src/*D.* *.log matrix* *~ *.pyc fluidity.* nohup.out *.node *.ele *.edge
