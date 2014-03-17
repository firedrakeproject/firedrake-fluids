import os
import pytest
import numpy
from firedrake import *

import les_smagorinsky

def mms_les_smagorinsky():
   errors = []
   for n in [4, 8, 16, 32]:
      mesh = UnitSquareMesh(n, n)
      smagorinsky_coefficient = 1.0*(1.0/n) # Remember to halve the Smagorinsky coefficient each time to compensate for the smaller filter width
      function_space = FunctionSpace(mesh, "CG", 1)
      u = [Function(function_space).interpolate(Expression('sin(x[0])')), Function(function_space).interpolate(Expression('cos(x[0])'))]
      smagorinsky_constant = Function(function_space).interpolate(Expression("1.0"))
      density = Function(function_space).interpolate(Expression("1.0"))
      exact_solution = Function(function_space).interpolate(Expression('%f * sqrt(2.0*cos(x[0])*cos(x[0]) + 0.5*sin(x[0])*sin(x[0]) + 0.5*sin(x[0])*sin(x[0]))' % smagorinsky_coefficient))
      visc = les_smagorinsky.eddy_viscosity(mesh, function_space, u, density, smagorinsky_coefficient, 1.0/n)
      print visc
      errors.append(sqrt(assemble(dot(visc - exact_solution, visc - exact_solution) * dx)))

   return errors

def test_mms_les_smagorinsky():
   errors = numpy.array(mms_les_smagorinsky())
   convergence_order = numpy.log2(errors[:-1] / errors[1:])
   print "Eddy viscosity convergence order:", convergence_order
   assert (numpy.array(convergence_order) > 0.8).all()
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
