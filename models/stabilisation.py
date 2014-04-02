# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

import fields_calculations

import numpy

def streamline_upwind(mesh, dimension, w, u, u_k, fs, h):

   scaling_factor = 0.5

   magnitude = fields_calculations.magnitude_vector(mesh, u_k, fs)

   u_nodes = magnitude.vector()
   near_zero = numpy.array([1.0e-9 for i in range(len(u_nodes))])
   u_nodes.set_local(numpy.maximum(u_nodes.array(), near_zero))

   k = 10.0 # Viscosity

   grid_pe = fields_calculations.grid_peclet_number(mesh, u_k, k, magnitude, fs, h)
   grid_pe_nodes = grid_pe.vector()

   values = numpy.array([1.0e-9 for i in range(len(grid_pe_nodes))])
   grid_pe_nodes.set_local(numpy.maximum(grid_pe_nodes.array(), values))

   k_bar = ( (1.0/tanh(grid_pe)) - (1.0/grid_pe) ) * h * magnitude

   F = 0
   for dim_i in range(dimension):
      for dim_j in range(dimension):
         F += 0.5*(k_bar/(magnitude**2))*inner(dot(grad(w[dim_i])[dim_j], u[dim_j]), dot(grad(u[dim_i])[dim_j], u[dim_j]))*dx

   return F

