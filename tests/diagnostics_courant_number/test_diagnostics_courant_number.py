import os
import pytest
import numpy
from firedrake import *

from diagnostics import *

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def diagnostics_courant_number():

   mesh = UnitSquareMesh(10, 10)
   fs = FunctionSpace(mesh, "CG", 1)
   velocity = [Function(fs).interpolate(Expression("x[0]")), Function(fs).interpolate(Expression("x[1]"))]
   dt = 0.5

   diagnostics = Diagnostics(mesh, fs)
   
   courant_number = diagnostics.courant_number(velocity, dt)
   courant_number_exact = Function(fs).interpolate(Expression("sqrt(pow(x[0], 2) + pow(x[1], 2))*0.5/0.1414213562"))
   norm = sqrt(assemble(dot(courant_number - courant_number_exact, courant_number - courant_number_exact) * dx))
   return norm

def test_diagnostics_courant_number(input):
   norm = diagnostics_courant_number()
   print "L2 error norm:", norm
   assert norm <= 9.0e-3
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
