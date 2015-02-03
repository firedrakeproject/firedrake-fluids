#    Copyright (C) 2013 Imperial College London.

#    This file is part of Firedrake-Fluids.
#
#    Firedrake-Fluids is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Firedrake-Fluids is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Firedrake-Fluids.  If not, see <http://www.gnu.org/licenses/>.

from firedrake import *

import firedrake_fluids.fields_calculations as fields_calculations

class Diagnostics:

   def __init__(self, mesh):
      self.mesh = mesh      
      return

   def courant_number(self, u, dt, function_space=None):
      r''' Compute the Courant number given by 
      
          .. math:: \frac{\|\mathbf{u}\|_2\Delta t}{\Delta x}
             
          where :math:`\mathbf{u}` is the velocity field, :math:`\Delta t` is the time-step, and :math:`\Delta x` is the element size.
          
      :param ufl.Function u: The velocity field.
      :param dt: The time-step size.
      :param ufl.FunctionSpace function_space: The function space in which the grid Reynolds number should be computed.
      :returns: A UFL Function representing the Courant number field.
      :rtype: ufl.Function
      '''
      
      if(function_space is None):
         function_space = FunctionSpace(self.mesh, "CG", 1)
         
      test = TestFunction(function_space)
      trial = TrialFunction(function_space)
      solution = Function(function_space)
      
      h = CellSize(self.mesh)
      magnitude = fields_calculations.magnitude_vector(u, function_space)
      
      a = inner(test, trial)*dx
      L = (test*(magnitude*dt)/h)*dx
      
      solve(a == L, solution, bcs=[])
      return solution
   
   def grid_reynolds_number(self, u, rho, mu=1, function_space=None):
      r''' Compute the grid Reynolds number given by 
      
          .. math:: \frac{\rho\|\mathbf{u}\|_2\Delta x}{\mu}
             
          where :math:`\rho` is the density, :math:`\mathbf{u}` is the velocity field, 
          :math:`\Delta x` is the element size, and :math:`\mu` is the (isotropic) viscosity.

      :param ufl.Function u: The velocity field.
      :param rho: The density, which can be a constant value or a ufl.Function.
      :param ufl.FunctionSpace function_space: The function space in which the grid Reynolds number should be computed.
      :param mu: The (isotropic) viscosity. By default, this is set to 1.0.
      :returns: A UFL Function representing the grid Reynolds number field.
      :rtype: ufl.Function
      '''

      if(function_space is None):
         function_space = FunctionSpace(self.mesh, "CG", 1)
         
      test = TestFunction(function_space)
      trial = TrialFunction(function_space)
      
      solution = Function(function_space)
      
      h = CellSize(self.mesh)
      magnitude = fields_calculations.magnitude_vector(u, function_space)
      
      a = inner(test, trial)*dx
      L = (test*(rho*magnitude*h)/mu)*dx
      
      solve(a == L, solution, bcs=[])
      return solution

   def divergence(self, u, function_space=None):
      r''' For a given vector field :math:`\mathbf{u}`, return its divergence (a scalar field):
      
           .. math:: \nabla\cdot\mathbf{u}
      
      :param ufl.Function u: A vector field.
      :param ufl.FunctionSpace function_space: The function space in which the grid Reynolds number should be computed.
      :returns: A UFL Function representing the divergence field.
      :rtype: ufl.Function
      '''

      if(function_space is None):
         function_space = FunctionSpace(self.mesh, "CG", 1)
         
      test = TestFunction(function_space)
      trial = TrialFunction(function_space)
      
      solution = Function(function_space)
      
      a = inner(test, trial)*dx
      L = test*div(u)*dx
      
      solve(a == L, solution, bcs=[])
      return solution

