# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

def magnitude_vector(mesh, u, function_space):
   """ Calculate the magnitude of a given vector 'u'. """

#   function_space = u[0].function_space() # Assumes all components of velocity live in the same function space
   w = TestFunction(function_space)
   magnitude = TrialFunction(function_space)
   solution = Function(function_space)

   a = w*magnitude*dx
   L = w*sqrt(dot(u, u))*dx
   solve(a == L, solution, bcs=[])

   return solution

def grid_peclet_number(mesh, diffusivity, magnitude, function_space, cellsize):
   """ Calculate the grid Peclet number, given by
       grid_pe = (|u|*cellsize)/(2*diffusivity) """

   w = TestFunction(function_space)
   grid_pe = TrialFunction(function_space)
   solution = Function(function_space)

   a = w*grid_pe*dx
   L = w*(magnitude*cellsize)/(2.0*diffusivity)*dx
   solve(a == L, solution, bcs=[])

   return solution

