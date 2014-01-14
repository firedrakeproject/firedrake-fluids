# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

import fields_calculations

def courant_number(mesh, u, dt):
   
   fs = u[0].function_space()
   w = TestFunction(fs)
   courant_number = TrialFunction(fs)
   solution = Function(fs)
   
   h = CellSize(mesh)
   magnitude = fields_calculations.magnitude_vector(u)
   
   a = inner(w, courant_number)*dx
   L = w*(magnitude*dt)/h*dx
   solve(a == L, solution, bcs=[])
   return solution
   
def divergence(vector):

   w = TestFunction(vector.function_space())
   divergence = TrialFunction(vector.function_space())
   solution = Function(vector.function_space())
   
   a = inner(w, divergence)*dx
   L = w*div(vector)*dx
   solve(a == L, solution, bcs=[])
   return solution

