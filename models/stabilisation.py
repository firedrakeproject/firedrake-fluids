#    Copyright (C) 2013 Imperial College London.

#    This file is part of Firedrake-Fluids.
#
#    Firedrake-Fluids is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Firedrake-Fluids is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Firedrake-Fluids.  If not, see <http://www.gnu.org/licenses/>.

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
      
