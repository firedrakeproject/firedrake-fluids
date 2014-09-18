import os
import pytest
import numpy
from firedrake import *

from firedrake_fluids.shallow_water import *

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def swe_wave_speed_1d():
   sw = ShallowWater(path=os.path.join(cwd, "swe_wave_speed_1d.swml"))
   sw.run()
   h_old = sw.solution_old.split()[-1]
   return h_old.vector().array()

def test_swe_wave_speed_1d(input):
   h_values = numpy.array(swe_wave_speed_1d())
   
   # Find the position of the wave front at the end of the simulation.
   dx = 1.0/1000.0
   wave_position = 0.0
   for i in range(0, len(h_values)):
      if(h_values[i] <= 0.5):
         wave_position = i*dx
         break
   
   # The wave should be travelling at sqrt(g*H) = 7 m/s.
   # Therefore, after 0.01 s, the wave should have travelled 0.07 m from the centre location at x = 0.5 m.
   print "Wave position = %f" % wave_position
   assert wave_position >= 0.43 - 10*dx and wave_position <= 0.43 + 10*dx

if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
