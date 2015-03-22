""" A test case described by Bermudez and Vazquez (1994). 
This involves tidal flow with a non-uniform bathymetry. 
The numerical results are compared with an asymptotic solution for
the velocity and free surface perturbation at t = 10,000 seconds. """

import os
import pytest
import numpy
import pylab
from firedrake import *

cwd = os.path.dirname(os.path.abspath(__file__))

def swe_tidal_flow_regular_bed():
   from firedrake_fluids.shallow_water import ShallowWater
   from firedrake_adjoint import project
   
   sw = ShallowWater(path=os.path.join(cwd, "swe_tidal_flow_regular_bed.swml"))
   solution = sw.run()
   
   mesh = sw.mesh
   
   # Get the coordinates of each node
   coordinates = mesh.coordinates
   
   # For projecting to the same space as the coordinate field
   vfs = VectorFunctionSpace(mesh, "CG", 1)
   fs = FunctionSpace(mesh, "CG", 1)
   
   u = solution.split()[0]
   h = solution.split()[1]
  
   u = project(u, vfs, annotate=False)
   h = project(h, fs, annotate=False)

   print coordinates.vector().array()
   print u.vector().array()
   print h.vector().array()
   
   # Use [0::2] to select the even-numbered components of the vector's array (i.e. the x-component only).
   return coordinates.vector().array()[0::2], u.vector().array()[0::2], h.vector().array()
   
def test_swe_tidal_flow_regular_bed():
   
   coordinates, ux, h = swe_tidal_flow_regular_bed()

   L = 14000.0
   t = 10000.0
   
   x = coordinates
   H = 50.5 - 40*x/L - 10*numpy.sin(numpy.pi*(4*x/L - 0.5))
   h_analytical = H + 4 - 4*numpy.sin(numpy.pi*(4*t/86400.0 + 0.5))
   ux_analytical = ((x - 14000.0)*numpy.pi/(5400.0*h_analytical))*numpy.cos(numpy.pi*(4*t/86400.0 + 0.5))
   
   plot = False
   if(plot):
      pylab.figure(0)
      pylab.plot(x, ux_analytical, 'g-', label="Analytical")
      pylab.plot(x[0:len(ux):5], ux[0:len(ux):5], 'b.', label="Numerical (Firedrake-Fluids)")
      pylab.legend()
      pylab.xlabel("x (m)")
      pylab.ylabel("Velocity (m/s)")
      pylab.savefig('swe_tidal_flow_regular_bed_velocity.png')
      
      pylab.figure(1)
      pylab.plot(x, h_analytical, 'g-', label="Analytical")
      pylab.plot(x[0:len(h):5], h[0:len(h):5] + H[0:len(h):5], 'b.', label="Numerical (Firedrake-Fluids)")
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
