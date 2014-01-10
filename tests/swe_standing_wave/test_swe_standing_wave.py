import os
import pytest
import numpy
from firedrake import *

import shallow_water

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def swe_standing_wave():
   sw = shallow_water.ShallowWater(path=os.path.join(cwd, "swe_standing_wave.swml"))
   sw.run()
   h_old = sw.h_old
   
   return h_old.vector().array()

def test_swe_standing_wave(input):
   h_values = numpy.array(swe_standing_wave())
   h_max = max(h_values)
   h_min = min(h_values)
   print "h_max: ", h_max
   print "h_min: ", h_min
   assert h_max <= 1.0
   assert h_min >= -1.0

if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
