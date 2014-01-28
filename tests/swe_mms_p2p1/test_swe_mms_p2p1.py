import os
import pytest
import numpy
from firedrake import *

import shallow_water

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def swe_mms_p2p1():
   configs = ["MMS_A", "MMS_B", "MMS_C"]
   ux_norms = []
   uy_norms = []
   h_norms = []
   
   for c in configs:
      sw = shallow_water.ShallowWater(path=os.path.join(cwd, c + ".swml"))
      sw.run()
      ux_old = sw.solution_old.split()[0]
      uy_old = sw.solution_old.split()[1]
      h_old = sw.solution_old.split()[-1]

      ux_exact = Function(sw.function_spaces["VelocityFunctionSpace"])
      ux_exact.interpolate(Expression("cos(x[0])"))
      uy_exact = Function(sw.function_spaces["VelocityFunctionSpace"])
      uy_exact.interpolate(Expression("0.0"))
      
      h_exact = Function(sw.function_spaces["FreeSurfaceFunctionSpace"])
      h_exact.interpolate(Expression("sin(x[0])*sin(x[1])"))
      
      h_norms.append(sqrt(assemble(dot(h_old - h_exact, h_old - h_exact) * dx)))
      ux_norms.append(sqrt(assemble(dot(ux_old - ux_exact, ux_old - ux_exact) * dx)))
      uy_norms.append(sqrt(assemble(dot(uy_old - uy_exact, uy_old - uy_exact) * dx)))

   return (h_norms, ux_norms, uy_norms)

def test_swe_mms_p2p1(input):
   h_norms, ux_norms, uy_norms = numpy.array(swe_mms_p2p1())
   h_convergence_order = numpy.log2(h_norms[:-1] / h_norms[1:])
   ux_convergence_order = numpy.log2(ux_norms[:-1] / ux_norms[1:])
   uy_convergence_order = numpy.log2(uy_norms[:-1] / uy_norms[1:])
   print "FreeSurfacePerturbation convergence order:", h_convergence_order
   print "Velocity (x-component) convergence order:", ux_convergence_order
   print "Velocity (y-component) convergence order:", uy_convergence_order
   assert (numpy.array(h_convergence_order) > 1.3).all()
   assert (numpy.array(ux_convergence_order) > 2.8).all()
   assert (numpy.array(uy_convergence_order) > 2.6).all()
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
