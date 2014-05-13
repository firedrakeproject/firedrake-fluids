""" A test case described by Bermudez and Vazquez (1994). 
This involves tidal flow with a non-uniform bathymetry. 
The numerical results are compared with an asymptotic solution for
the velocity and free surface perturbation at t = 10,000 seconds. """

import os
import pytest
import numpy
import pylab
from firedrake import *

import shallow_water

cwd = os.path.dirname(os.path.abspath(__file__))

def swe_tidal_flow_regular_bed():
   sw = shallow_water.ShallowWater(path=os.path.join(cwd, "swe_tidal_flow_regular_bed.swml"))
   sw.run()
   
   mesh = sw.mesh
   
   # Get the coordinates of each node
   coordinates = mesh.coordinates
   
   # For projecting to the same space as the coordinate field
   fs = FunctionSpace(sw.mesh, "CG", 1)
   
   u_old = sw.solution_old.split()[0]
   h_old = sw.solution_old.split()[1]
   
   u_old = project(u_old[0], fs)
   h_old = project(h_old, fs)
   
   return coordinates.vector().array(), u_old.vector().array(), h_old.vector().array()
   
def test_swe_tidal_flow_regular_bed():
   
   coordinates, ux, h = swe_tidal_flow_regular_bed()

   L = 14000.0
   t = 10000.0
   
   x = coordinates[:, 0]
   H = 50.5 - 40*x/L - 10*numpy.sin(numpy.pi*(4*x/L - 0.5))
   h_analytical = H + 4 - 4*numpy.sin(numpy.pi*(4*t/86400.0 + 0.5))
   ux_analytical = ((x - 14000.0)*numpy.pi/(5400.0*h_analytical))*numpy.cos(numpy.pi*(4*t/86400.0 + 0.5))
   
   plot = False
   if(plot):
      pylab.figure(0)
      pylab.plot(x, ux_analytical, 'g-', label="Analytical")
      pylab.plot(x[0:len(ux):5], ux[0:len(ux):5], 'b.', label="Numerical (Firedrake)")
      pylab.legend()
      pylab.xlabel("x (m)")
      pylab.ylabel("Velocity (m/s)")
      pylab.savefig('swe_tidal_flow_regular_bed_velocity.png')
      
      pylab.figure(1)
      pylab.plot(x, h_analytical, 'g-', label="Analytical")
      pylab.plot(x[0:len(h):5], h[0:len(h):5] + H[0:len(h):5], 'b.', label="Numerical (Firedrake)")
      pylab.legend()
      pylab.xlabel("x (m)")
      pylab.ylabel("Free surface (m)")
      pylab.axis([0, 14000, 0, 70])
      pylab.savefig('swe_tidal_flow_regular_bed_freesurface.png')

   # Check that the maximum error between the numerical and analytical solution is small.
   ux_difference = abs(ux_analytical - ux)
   assert(max(ux_difference) < 5.0e-3)
   
   height = numpy.zeros(len(h))
   for i in range(len(h)):
      height[i] = h[i] + H[i]
   h_difference = abs(h_analytical - height)
   assert(max(h_difference) < 5.0e-3)

if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
