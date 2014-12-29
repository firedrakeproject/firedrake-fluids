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

   def __init__(self, mesh, function_space):
      self.mesh = mesh
      self.function_space = function_space
      
      self.test = TestFunction(self.function_space)
      self.trial = TrialFunction(self.function_space)
      
      return

   def courant_number(self, u, dt):
      r''' Compute the Courant number given by 
      
          .. math:: \frac{\|\mathbf{u}\|_2\Delta t}{\Delta x}
             
          where :math:`\mathbf{u}` is the velocity field, :math:`\Delta t` is the time-step, and :math:`\Delta x` is the element size.
          
      :param ufl.Function u: The velocity field.
      :param dt: The time-step size.
      :returns: A UFL Function representing the Courant number field.
      :rtype: ufl.Function
      '''
      solution = Function(self.function_space)
      
      h = CellSize(self.mesh)
      magnitude = fields_calculations.magnitude_vector(u, self.function_space)
      
      a = inner(self.test, self.trial)*dx
      L = (self.test*(magnitude*dt)/h)*dx
      
      solve(a == L, solution, bcs=[])
      return solution
   
   def grid_reynolds_number(self, rho, u, mu):
      r''' Compute the grid Reynolds number given by 
      
          .. math:: \frac{\rho\|\mathbf{u}\|_2\Delta x}{\mu}
             
          where :math:`\rho` is the density, :math:`\mathbf{u}` is the velocity field, 
          :math:`\Delta x` is the element size, and :math:`\mu` is the (isotropic) viscosity.
         
      :param rho: The density, which can be a constant value or a ufl.Function.
      :param ufl.Function u: The velocity field.
      :param mu: The (isotropic) viscosity.
      :returns: A UFL Function representing the grid Reynolds number field.
      :rtype: ufl.Function
      '''
      solution = Function(self.function_space)
      
      h = CellSize(self.mesh)
      magnitude = fields_calculations.magnitude_vector(u, self.function_space)
      
      a = inner(self.test, self.trial)*dx
      L = (self.test*(rho*magnitude*h)/mu)*dx
      
      solve(a == L, solution, bcs=[])
      return solution

   def divergence(self, u):
      r''' For a given vector field :math:`\mathbf{u}`, return its divergence (a scalar field):
      
           .. math:: \nabla\cdot\mathbf{u}
      
      :param ufl.Function u: A vector field.
      :returns: A UFL Function representing the divergence field.
      :rtype: ufl.Function
      '''
      
      solution = Function(self.function_space)
      
      a = inner(self.test, self.trial)*dx
      L = self.test*div(u)*dx
      
      solve(a == L, solution, bcs=[])
      return solution

