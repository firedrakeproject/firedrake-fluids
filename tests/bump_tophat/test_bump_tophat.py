import os
import pytest
import numpy
from firedrake import *
from firedrake_fluids.shallow_water import *

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)
   
def bump_tophat():

   sw = ShallowWater(path=os.path.join(cwd, "bump.swml"))
   sw.run()
   bump = assemble(sw.array.turbine_drag*dx)

   sw = ShallowWater(path=os.path.join(cwd, "tophat.swml"))
   sw.run()
   tophat = assemble(sw.array.turbine_drag*dx)
   
   return bump, tophat

def test_bump_tophat(input):
   bump, tophat = bump_tophat()
   assert abs(tophat - 1250) < 50.0
   assert abs(bump - 1250) < 50.0
   assert abs(tophat - bump) < 20.0
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
