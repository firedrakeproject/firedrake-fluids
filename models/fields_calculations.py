# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

def magnitude_vector(u):

   function_space = u[0].function_space() # Assumes all components of velocity live in the same function space
   w = TestFunction(function_space)
   magnitude = TrialFunction(function_space)
   solution = Function(function_space)
   
   a = w*magnitude*dx
   L = w*sqrt(sum([dot(u[i], u[i]) for i in range(0, len(u))]))*dx
   solve(a == L, solution, bcs=[])
   
   return solution


