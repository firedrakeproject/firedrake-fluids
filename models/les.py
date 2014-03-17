# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *
from random import random

def element_volume(v):
   # FIXME: There is currently no way to get the cell volume in Firedrake.
   return v.cell().volume()
   
def strain_rate_tensor(u):
   dimension = len(u)
   S = [[0 for j in range(dimension)] for i in range(dimension)]
   for dim_i in range(dimension):
      for dim_j in range(dimension):
         S[dim_i][dim_j] = -0.5*(grad(u[dim_i])[dim_j] + grad(u[dim_j])[dim_i])
   return S

def eddy_viscosity(mesh, function_space, u, density, smagorinsky_coefficient, filter_width):

   dimension = len(u)
   w = TestFunction(function_space)
   eddy_viscosity = TrialFunction(function_space)

   #if(dimension == 2):
   #   filter_width = element_volume(eddy_viscosity)**(1.0/2.0)
   #elif(dimension == 3):
   #   filter_width = element_volume(eddy_viscosity)**(1.0/3.0)
   #else:
   #   print "Dimension == 1"
      
   S = strain_rate_tensor(u)
   second_invariant = 0.0
   for i in range(0, dimension):
      for j in range(0, dimension):
         second_invariant += 2.0*(S[i][j]**2)
   second_invariant = sqrt(second_invariant)
   rhs = density*(smagorinsky_coefficient*filter_width)**2*second_invariant

   a = inner(w, eddy_viscosity)*dx
   L = inner(w, rhs)*dx

   solution = Function(function_space)
   solve(a == L, solution, bcs=[])

   return solution

