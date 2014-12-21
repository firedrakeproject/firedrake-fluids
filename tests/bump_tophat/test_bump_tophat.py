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
   expected_tophat_integral = 20.0*5.0*12.0
   expected_bump_integral = 1.45661*8.5*10*10 # The 1.45661 is the area under a bump function with x = (-1, 1) and y = (-1, 1). So we need to multiply this by 10 m * 10 m.
   print "Expected tophat integral: %.2f" % expected_tophat_integral
   print "Computed tophat integral: %.2f" % tophat
   print "Expected bump integral: %.2f" % expected_bump_integral
   print "Computed bump integral: %.2f" % bump

   assert abs(bump - expected_bump_integral) <= 0.5
   assert abs(tophat - expected_tophat_integral) <= 55.0
   assert abs(bump - expected_tophat_integral) <= 55.0
   assert abs(tophat - bump) < 20.0
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
