import os
import pytest
import numpy
from firedrake import *

from diagnostics import *

cwd = os.path.dirname(os.path.abspath(__file__))

@pytest.fixture(scope='session')
def input():
   os.system("make -C " + cwd)

def diagnostics_grid_reynolds_number():

   mesh = UnitSquareMesh(10, 10)
   fs = FunctionSpace(mesh, "CG", 1)
   velocity = [Function(fs).interpolate(Expression("x[0]")), Function(fs).interpolate(Expression("x[1]"))]
   density = Function(fs).interpolate(Expression("2.0"))
   mu = Function(fs).interpolate(Expression("0.7"))

   diagnostics = Diagnostics(mesh, fs)

   grid_re_number = diagnostics.grid_reynolds_number(density, velocity, mu)
   grid_re_number_exact = Function(fs).interpolate(Expression("2.0*sqrt(pow(x[0], 2) + pow(x[1], 2))*0.1414213562/0.7"))
   norm = sqrt(assemble(dot(grid_re_number - grid_re_number_exact, grid_re_number - grid_re_number_exact) * dx))
   return norm

def test_diagnostics_grid_reynolds_number(input):
   norm = diagnostics_grid_reynolds_number()
   print "L2 error norm:", norm
   assert norm <= 5.0e-4
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
