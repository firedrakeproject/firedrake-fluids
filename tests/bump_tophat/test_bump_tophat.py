import os
import pytest
import numpy
from firedrake import *
from firedrake_fluids.shallow_water import *

cwd = os.path.dirname(os.path.abspath(__file__))

def bump_tophat():

   sw = ShallowWater(path=os.path.join(cwd, "bump.swml"))
   sw.run()
   bump = assemble(sw.array.turbine_drag*dx)

   sw = ShallowWater(path=os.path.join(cwd, "tophat.swml"))
   sw.run()
   tophat = assemble(sw.array.turbine_drag*dx)
   
   return bump, tophat

def test_les_smagorinsky_eddy_viscosity():
   bump, tophat = bump_tophat()
   assert abs(tophat - bump) < 100.0
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
