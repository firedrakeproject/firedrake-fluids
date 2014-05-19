# Copyright 2013 Imperial College London. All rights reserved.

from firedrake import *

import fields_calculations

import numpy

class Stabilisation:

   def __init__(self, mesh, function_space, cellsize):
      self.mesh = mesh
      self.function_space = function_space
      self.cellsize = cellsize

   def streamline_upwind(self, w, u, magnitude, grid_pe):

      dimension = len(u)

      scaling_factor = 0.5
      k_bar = self.k_bar(magnitude, grid_pe)

      F = scaling_factor*(k_bar/(magnitude**2))*inner(dot(grad(w), u), dot(grad(u), u))*dx

      return F

   def k_bar(self, magnitude, grid_pe):

      return ( (1.0/tanh(grid_pe)) - (1.0/grid_pe) ) * self.cellsize * magnitude
      
