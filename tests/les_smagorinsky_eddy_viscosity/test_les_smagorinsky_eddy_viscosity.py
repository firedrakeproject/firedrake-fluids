import os
import pytest
import numpy
from firedrake import *

from les import *

def les_smagorinsky_eddy_viscosity():
   errors = []
   for n in [2, 4, 8, 16, 32]:
      mesh = UnitSquareMesh(n, n)
      smagorinsky_coefficient = 2.0
      filter_width = sqrt((1.0/n)**2 + (1.0/n)**2) # FIXME: Use CellSize instead.
      
      fs_exact = FunctionSpace(mesh, "CG", 3)
      fs = FunctionSpace(mesh, "CG", 1)
      
      u = [Function(fs).interpolate(Expression('sin(x[0])')), Function(fs).interpolate(Expression('cos(x[0])'))]
      density = Function(fs).interpolate(Expression("1.0"))
      
      exact_solution = project(Expression('pow(%f, 2) * sqrt(2.0*cos(x[0])*cos(x[0]) + 0.5*sin(x[0])*sin(x[0]) + 0.5*sin(x[0])*sin(x[0]))' % smagorinsky_coefficient), fs_exact)

      les = LES(mesh, fs)
      # Since the eddy viscosity depends on the filter_width, we need to provide smagorinsky_coefficient/filter_width here 
      # for the convergence test because we want to compare against the same exact solution throughout.
      eddy_viscosity = les.eddy_viscosity(u, density, (smagorinsky_coefficient/filter_width), filter_width)
      print eddy_viscosity.vector()[:]
      errors.append(sqrt(assemble(dot(eddy_viscosity - exact_solution, eddy_viscosity - exact_solution) * dx)))

   return errors

def test_les_smagorinsky_eddy_viscosity():
   errors = numpy.array(les_smagorinsky_eddy_viscosity())
   assert (numpy.array([numpy.log2(errors[i]/errors[i+1]) for i in range(len(errors)-1)]) > 1.45).all()
   
if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
