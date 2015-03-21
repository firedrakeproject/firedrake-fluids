import os
import pytest
import numpy
from firedrake import *

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def swe_dam_break_2d():
   from firedrake_fluids.shallow_water import ShallowWater
   sw = ShallowWater(path=os.path.join(cwd, "swe_dam_break_2d.swml"))
   solution = sw.run()
   
   h = solution.split()[1]
   u = solution.split()[0]
   
   return h.vector().array(), u.vector().array()

def test_swe_dam_break_2d(input):
   h, u = swe_dam_break_2d()
   
   # x-component of Velocity
   assert(max(u[0::2]) <= 9.0)   
   assert(max(u[0::2]) >= -2.0)
   # y-component of Velocity
   assert(max(u[1::2]) <= 6.0)   
   assert(max(u[1::2]) >= -6.0)
   # Free surface
   assert(max(h) <= 5.0)   
   assert(max(h) >= -4.0)

if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
