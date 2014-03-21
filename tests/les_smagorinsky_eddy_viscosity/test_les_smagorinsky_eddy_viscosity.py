import os
import pytest
import numpy
from firedrake import *

from les import *

def les_smagorinsky_eddy_viscosity():
   errors = []
   for n in [4, 8, 16, 32]:
      mesh = UnitSquareMesh(n, n)
      smagorinsky_coefficient = 1.0
      filter_width = sqrt((1.0/n)**2 + (1.0/n)**2) # FIXME: Use CellSize instead.
      function_space = FunctionSpace(mesh, "CG", 1)
      u = [Function(function_space).interpolate(Expression('sin(x[0])')), Function(function_space).interpolate(Expression('cos(x[0])'))]
      density = Function(function_space).interpolate(Expression("1.0"))
      exact_solution = Function(function_space).interpolate(Expression('pow(%f*%f, 2) * sqrt(2.0*cos(x[0])*cos(x[0]) + 0.5*sin(x[0])*sin(x[0]) + 0.5*sin(x[0])*sin(x[0]))' % (smagorinsky_coefficient, filter_width)))

      les = LES(mesh, function_space)
      visc = les.eddy_viscosity(u, density, smagorinsky_coefficient, 1.0/n)
      print visc.vector()[:]
      errors.append(sqrt(assemble(dot(visc - exact_solution, visc - exact_solution) * dx)))

   return errors

def test_les_smagorinsky_eddy_viscosity():
   errors = numpy.array(les_smagorinsky_eddy_viscosity())
   convergence_order = numpy.log2(errors[:-1] / errors[1:])
   print "Eddy viscosity convergence order:", convergence_order
   assert (numpy.array(convergence_order) > 1.9).all()
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
