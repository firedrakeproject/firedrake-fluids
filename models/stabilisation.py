# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

import fields_calculations

def streamline_upwind(mesh, dimension, w, u, u_k):

   h = CellSize(mesh)
   scaling_factor = 0.5
   
   magnitude = fields_calculations.magnitude_vector(u_k)

   u_nodes = magnitude.vector()
   for i in range(0, len(u_nodes)):
      if(u_nodes[i] < 1.0e-9):
         u_nodes[i] = 1.0e-9

   #k = 10.0 # Viscosity
   #grid_pe = (magnitude*h)/(2.0*k)
   #k_bar = ( (1.0/tanh(grid_pe)) - (1.0/grid_pe) ) * k * grid_pe

   k_bar = scaling_factor*h*magnitude

   F = 0
   for dim in range(dimension):
      F += (k_bar/(magnitude**2))*inner(u_k[dim], u_k[dim])*inner(grad(u[dim]), grad(w[dim]))*dx

   return F

