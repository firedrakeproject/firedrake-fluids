import os
import pytest
import numpy
from firedrake import *

cwd = os.path.dirname(os.path.abspath(__file__))

def swe_steady_state():
   from firedrake_fluids.shallow_water import ShallowWater
   sw = ShallowWater(path=os.path.join(cwd, "swe_steady_state.swml"))
   solution = sw.run()
   u_old = solution.split()[0]
   h_old = solution.split()[1]
   
   return u_old.vector().array(), h_old.vector().array()
   
def test_swe_steady_state():
   
   u_values, h_values = numpy.array(swe_steady_state())
   ux_values = u_values[0::2]
   uy_values = u_values[1::2]
   
   ux_max = max(ux_values)
   ux_min = min(ux_values)
   uy_max = max(uy_values)
   uy_min = min(uy_values)
   h_max = max(h_values)
   h_min = min(h_values)
   
   print "ux_max: ", ux_max
   print "ux_min: ", ux_min
   print "uy_max: ", uy_max
   print "uy_min: ", uy_min
   print "h_max: ", h_max
   print "h_min: ", h_min
   
   # A steady-state Velocity (in the x direction) of 2.0 is desired.
   assert ux_max <= 2.01
   assert ux_min >= 1.99
   
   # The y component of the Velocity field should be close to zero in the channel.
   assert uy_max <= 1e-4
   assert uy_min >= -1e-4
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
