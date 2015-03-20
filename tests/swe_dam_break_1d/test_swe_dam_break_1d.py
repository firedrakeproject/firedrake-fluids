import os
import pytest
import numpy
from firedrake import *

from firedrake_fluids.shallow_water import *

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def swe_dam_break_1d():
   sw = ShallowWater(path=os.path.join(cwd, "swe_dam_break_1d.swml"))
   solution = sw.run()
   h_old = solution.split()[1]
   u_old = solution.split()[0]
   
   return h_old.vector().array(), u_old.vector().array()

def test_swe_dam_break_1d(input):
   h, u = swe_dam_break_1d()
   
   assert(max(u) <= 3.0 + 1.0e-2)   
   assert(max(u) >= 0.0 - 1.0e-2)
   assert(max(h) <= 5.0 + 1.0e-2)   
   assert(max(h) >= 0.0 - 1.0e-2)

if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
