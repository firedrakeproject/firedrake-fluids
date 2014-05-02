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
   
   fs = FunctionSpace(sw.mesh, "CG", 1)
   vfs = VectorFunctionSpace(sw.mesh, "CG", 2)
   
   h_old = project(sw.solution_old.split()[1], fs)
   u_old = project(sw.solution_old.split()[0], vfs)   
   
   return h_old.vector().array(), u_old.vector().array()

def test_swe_dam_break_2d(input):
   h, u = swe_dam_break_2d()
   
   assert(max(u[:,0]) <= 9.0)   
   assert(max(u[:,0]) >= -2.0)
   assert(max(u[:,1]) <= 6.0)   
   assert(max(u[:,1]) >= -6.0)
   assert(max(h) <= 5.0)   
   assert(max(h) >= -4.0)

if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
