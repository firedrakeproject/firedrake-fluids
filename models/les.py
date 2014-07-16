# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *
from random import random

class LES:

   def __init__(self, mesh, function_space):
      self.mesh = mesh
      self.function_space = function_space
      return
      
   def strain_rate_tensor(self, u):
      S = 0.5*(grad(u) + grad(u).T)
      return S

   def eddy_viscosity(self, u, density, smagorinsky_coefficient):

      dimension = len(u)
      w = TestFunction(self.function_space)
      eddy_viscosity = TrialFunction(self.function_space)

      filter_width = CellVolume(self.mesh)**(1.0/dimension)

      S = self.strain_rate_tensor(u)
      second_invariant = 0.0
      for i in range(0, dimension):
         for j in range(0, dimension):
            second_invariant += 2.0*(S[i,j]**2)
            
      second_invariant = sqrt(second_invariant)
      rhs = density*(smagorinsky_coefficient*filter_width)**2*second_invariant

      lhs = inner(w, eddy_viscosity)*dx
      rhs = inner(w, rhs)*dx
      
      return lhs, rhs

