import os
import pytest
import numpy
from firedrake import *

from firedrake_fluids.shallow_water import *

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def swe_mms_p2p1_quadratic_drag_divergence_free():
   configs = ["MMS_A", "MMS_B", "MMS_C"]
   ux_norms = []
   uy_norms = []
   h_norms = []
   
   for c in configs:
      sw = ShallowWater(path=os.path.join(cwd, c + ".swml"))
      solution = sw.run()
      ux_old = solution.split()[0][0]
      uy_old = solution.split()[0][1]
      h_old = solution.split()[1]
      
      fs_exact = FunctionSpace(sw.mesh, "CG", 3)
      
      ux_exact = project(Expression("cos(x[1])*sin(x[0])"), fs_exact)
      uy_exact = project(Expression("-cos(x[0])*sin(x[1])"), fs_exact)
      h_exact = project(Expression("sin(x[0])*sin(x[1])"), fs_exact)
      
      h_norms.append(sqrt(assemble(dot(h_old - h_exact, h_old - h_exact) * dx)))
      ux_norms.append(sqrt(assemble(dot(ux_old - ux_exact, ux_old - ux_exact) * dx)))
      uy_norms.append(sqrt(assemble(dot(uy_old - uy_exact, uy_old - uy_exact) * dx)))

   return (h_norms, ux_norms, uy_norms)

def test_swe_mms_p2p1_quadratic_drag_divergence_free(input):
   h_norms, ux_norms, uy_norms = numpy.array(swe_mms_p2p1_quadratic_drag_divergence_free())
   print h_norms
   print ux_norms
   print uy_norms
   print numpy.array([numpy.log2(h_norms[i]/h_norms[i+1]) for i in range(len(h_norms)-1)])
   print numpy.array([numpy.log2(ux_norms[i]/ux_norms[i+1]) for i in range(len(ux_norms)-1)])
   print numpy.array([numpy.log2(uy_norms[i]/uy_norms[i+1]) for i in range(len(uy_norms)-1)])

   assert (numpy.array([numpy.log2(h_norms[i]/h_norms[i+1]) for i in range(len(h_norms)-1)]) > 2.0).all()
   assert (numpy.array([numpy.log2(ux_norms[i]/ux_norms[i+1]) for i in range(len(ux_norms)-1)]) > 2.65).all()
   assert (numpy.array([numpy.log2(uy_norms[i]/uy_norms[i+1]) for i in range(len(uy_norms)-1)]) > 2.7).all()
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
