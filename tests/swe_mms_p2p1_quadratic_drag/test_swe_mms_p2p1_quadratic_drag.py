import os
import pytest
import numpy
from firedrake import *

import shallow_water

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def swe_mms_p2p1_quadratic_drag():
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
      
      fs_exact = FunctionSpace(sw.mesh, "CG", 3)
      
      ux_exact = project(Expression("pow(cos(x[1]), 3)*sin(x[0])"), fs_exact)
      uy_exact = project(Expression("pow(sin(x[0]), 2)*sin(x[1])"), fs_exact)
      h_exact = project(Expression("sin(x[0])*sin(x[1])"), fs_exact)
      
      h_norms.append(sqrt(assemble(dot(h_old - h_exact, h_old - h_exact) * dx)))
      ux_norms.append(sqrt(assemble(dot(ux_old - ux_exact, ux_old - ux_exact) * dx)))
      uy_norms.append(sqrt(assemble(dot(uy_old - uy_exact, uy_old - uy_exact) * dx)))

   return (h_norms, ux_norms, uy_norms)

def test_swe_mms_p2p1_quadratic_drag(input):
   h_norms, ux_norms, uy_norms = numpy.array(swe_mms_p2p1_quadratic_drag())
   assert (numpy.array([numpy.log2(h_norms[i]/h_norms[i+1]) for i in range(len(h_norms)-1)]) > 1.45).all()
   assert (numpy.array([numpy.log2(ux_norms[i]/ux_norms[i+1]) for i in range(len(ux_norms)-1)]) > 2.2).all()
   assert (numpy.array([numpy.log2(uy_norms[i]/uy_norms[i+1]) for i in range(len(uy_norms)-1)]) > 1.4).all()
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
