# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

import fields_calculations

import numpy

class Stabilisation:

   def __init__(self, mesh, function_space, cellsize):
      self.mesh = mesh
      self.function_space = function_space
      self.cellsize = cellsize

   def streamline_upwind(self, w, u, u_old, diffusivity):

      dimension = len(u)
      
      magnitude = fields_calculations.magnitude_vector(self.mesh, u_old, self.function_space)
      
      # Bound the values for the magnitude below by 1.0e-9 for numerical stability reasons. 
      u_nodes = magnitude.vector()
      near_zero = numpy.array([1.0e-9 for i in range(len(u_nodes))])
      u_nodes.set_local(numpy.maximum(u_nodes.array(), near_zero))
      
      scaling_factor = 0.5
      k_bar = self.k_bar(magnitude, diffusivity)

      F = 0
      for dim_i in range(dimension):
         for dim_j in range(dimension):
            F += scaling_factor*(k_bar/(magnitude**2))*inner(dot(grad(w[dim_i])[dim_j], u[dim_j]), dot(grad(u[dim_i])[dim_j], u[dim_j]))*dx

      return F

   def k_bar(self, magnitude, diffusivity):

      grid_pe = fields_calculations.grid_peclet_number(self.mesh, diffusivity, magnitude, self.function_space, self.cellsize)
      
      # Bound the values for grid_pe below by 1.0e-9 for numerical stability reasons. 
      grid_pe_nodes = grid_pe.vector()
      values = numpy.array([1.0e-9 for i in range(len(grid_pe_nodes))])
      grid_pe_nodes.set_local(numpy.maximum(grid_pe_nodes.array(), values))

      k_bar = ( (1.0/tanh(grid_pe)) - (1.0/grid_pe) ) * self.cellsize * magnitude
      
      return k_bar
      
