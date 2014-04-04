import os
import pytest
import numpy
from firedrake import *

from stabilisation import *

def su_stabilisation_k_bar():
   errors = []
   for n in [4, 8, 16, 32]:
      mesh = UnitSquareMesh(n, n)
      function_space = FunctionSpace(mesh, "CG", 1)
      
      viscosity = 0.1
      magnitude = Function(function_space).interpolate(Expression("sqrt(pow(x[0], 2) + pow(x[1], 2))"))

      cellsize_exact = sqrt((1.0/n)**2 + (1.0/n)**2)
      k_bar_exact = Function(function_space).interpolate(Expression("sqrt(pow(x[0], 2) + pow(x[1], 2))*(-2.00*0.1/(sqrt(pow(x[0], 2) + pow(x[1], 2))) + %f*1/tanh(0.500*sqrt(pow(x[0], 2) + pow(x[1], 2))*%f/0.1))" % (cellsize_exact, cellsize_exact)))

      stabilisation = Stabilisation(mesh, function_space, CellSize(mesh))
      k_bar = stabilisation.k_bar(magnitude, viscosity)
      
      test = TestFunction(function_space)
      trial = TrialFunction(function_space)
      solution = Function(function_space)
      solve(inner(test, trial)*dx == test*k_bar*dx, solution, bcs=[])
      
      errors.append(sqrt(assemble(dot(k_bar - k_bar_exact, k_bar - k_bar_exact) * dx)))

   return errors

def test_su_stabilisation_k_bar():
   """ Tests the correctness of the k_bar diffusivity term in the SU stabilisation term. """
   errors = numpy.array(su_stabilisation_k_bar())
   convergence_order = numpy.log2(errors[:-1] / errors[1:])
   print "Convergence order for the k_bar term in SU stabilisation:", convergence_order
   assert (numpy.array(convergence_order) > 1.8).all()
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
