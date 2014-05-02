import os
import pytest
import numpy
from firedrake import *

import shallow_water

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def swe_dam_break_2d():
   sw = shallow_water.ShallowWater(path=os.path.join(cwd, "swe_dam_break_2d.swml"))
   sw.run()
   h_old = sw.solution_old.split()[1]
   ux_old = sw.solution_old.split()[0][0]
   uy_old = sw.solution_old.split()[0][1]
   
   return h_old.vector().array(), ux_old.vector().array(), uy_old.vector().array()

def test_swe_dam_break_2d(input):
   h, ux, uy = swe_dam_break_2d()
   
   assert(max(ux) <= 9.0)   
   assert(max(ux) >= -2.0)
   assert(max(uy) <= 6.0)   
   assert(max(uy) >= -6.0)
   assert(max(h) <= 5.0)   
   assert(max(h) >= -4.0)

if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
