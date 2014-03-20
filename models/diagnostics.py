# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

import fields_calculations

class Diagnostics:

   def __init__(self, mesh, function_space):
      self.mesh = mesh
      self.function_space = function_space
      
      self.test = TestFunction(self.function_space)
      self.trial = TrialFunction(self.function_space)
      
      return

   def courant_number(self, u, dt):
      ''' Compute the Courant number given by 
             (u*dt)/dx
          where u is the velocity, dt is the time-step, and dx is the characteristic element length. '''
      solution = Function(self.function_space)
      
      h = CellSize(self.mesh)
      magnitude = fields_calculations.magnitude_vector(self.mesh, u)
      
      a = inner(self.test, self.trial)*dx
      L = (self.test*(magnitude*dt)/h)*dx
      
      solve(a == L, solution, bcs=[])
      return solution
   
   def grid_reynolds_number(self, rho, u, mu):
      ''' Compute the grid Reynolds number given by 
             (rho*u*dx)/mu
          where rho is the density, u is the velocity, 
          dx is the characteristic element length, and mu is the viscosity. '''
          
      solution = Function(self.function_space)
      
      h = CellSize(self.mesh)
      magnitude = fields_calculations.magnitude_vector(self.mesh, u)
      
      a = inner(self.test, self.trial)*dx
      L = (self.test*(rho*magnitude*h)/mu)*dx
      
      solve(a == L, solution, bcs=[])
      return solution

   def divergence(self, u):
      ''' Compute the divergence of a vector field u. Returns a scalar quantity. '''
      
      solution = Function(self.function_space)
      
      a = inner(self.test, self.trial)*dx
      L = self.test*div(u)*dx
      
      solve(a == L, solution, bcs=[])
      return solution

