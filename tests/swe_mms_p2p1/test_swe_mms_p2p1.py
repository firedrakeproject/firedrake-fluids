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
   norms = []

   for c in configs:
      sw = shallow_water.ShallowWater(path=os.path.join(cwd, c + ".swml"))
      sw.run()
      h_old = sw.solution_old.split()[-1]

      f = Function(sw.function_spaces["FreeSurfaceFunctionSpace"])
      f.interpolate(Expression("sin(x[0])*sin(x[1])"))

      norms.append(sqrt(assemble(dot(h_old - f, h_old - f) * dx)))

   return norms

def test_swe_mms_p2p1(input):
   norms = numpy.array(swe_mms_p2p1())
   convergence_order = numpy.log2(norms[:-1] / norms[1:])
   print "Convergence order:", convergence_order
   assert (numpy.array(convergence_order) > 1.3).all()

if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
