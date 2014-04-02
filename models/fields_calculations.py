# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

def magnitude_vector(mesh, u, function_space):

#   function_space = u[0].function_space() # Assumes all components of velocity live in the same function space
   w = TestFunction(function_space)
   magnitude = TrialFunction(function_space)
   solution = Function(function_space)

   a = w*magnitude*dx
   L = w*sqrt(sum([dot(u[i], u[i]) for i in range(0, len(u))]))*dx
   solve(a == L, solution, bcs=[])

   return solution

def grid_peclet_number(mesh, u, k, magnitude, function_space, h):

   w = TestFunction(function_space)
   grid_pe = TrialFunction(function_space)
   solution = Function(function_space)

   a = w*grid_pe*dx
   L = w*(magnitude*h)/(2.0*k)*dx
   solve(a == L, solution, bcs=[])

   return solution

