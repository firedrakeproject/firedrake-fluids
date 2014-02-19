# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

import fields_calculations

def streamline_upwind(mesh, dimension, w, u, u_k):

   h = CellSize(mesh)
   scaling_factor = 0.5
   
   magnitude = fields_calculations.magnitude_vector(mesh, u_k)

   u_nodes = magnitude.vector()
   for i in range(0, len(u_nodes)):
      if(u_nodes[i] < 1.0e-16):
         u_nodes[i] = 1.0e-16

   k = 10.0 # Viscosity

   grid_pe = fields_calculations.grid_peclet_number(mesh, u_k, k, magnitude)   
   grid_pe_nodes = grid_pe.vector()
   for i in range(0, len(grid_pe_nodes)):
      if(grid_pe_nodes[i] < 1.0e-16):
         grid_pe_nodes[i] = 1.0e-16
   for i in range(0, len(grid_pe_nodes)):
      if(grid_pe_nodes[i] < -1.0):
         grid_pe_nodes[i] = -1.0 - 1.0/grid_pe_nodes[i]
      elif(grid_pe_nodes[i] >= -1.0 and grid_pe_nodes[i] <= 1.0):
         grid_pe_nodes[i] = 0.0
      else:
         grid_pe_nodes[i] = 1.0 + 1.0/grid_pe_nodes[i]

   #grid_pe = (magnitude*h)/(2.0*k)

   #k_bar = grid_pe*h

   #grid_pe = (magnitude*h)/(2.0*k)
#   k_bar = ( (1.0/tanh(grid_pe)) - (1.0/grid_pe) ) * h
   k_bar = grid_pe*h*magnitude

   F = 0
   for dim_i in range(dimension):
      for dim_j in range(dimension):
         F += 0.5*(k_bar/(magnitude**2))*inner(dot(grad(w[dim_i])[dim_j], u[dim_j]), dot(grad(u[dim_i])[dim_j], u[dim_j]))*dx

   #F = 0
   #for dim in range(dimension):
   #   F += 0.5*(k_bar/(magnitude**2))*inner(u[dim], u[dim])*inner(grad(u[dim]), grad(w[dim]))*dx

   return F

