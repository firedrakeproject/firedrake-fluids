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

class Stabilisation:
   """ Module for advection stabilisation terms. """

   def __init__(self, mesh, function_space, cellsize):
      self.mesh = mesh
      self.function_space = function_space
      self.cellsize = cellsize

   def streamline_upwind(self, w, u, magnitude, grid_pe):
      r""" Construct the form for the streamline-upwind stabilisation term, given by equation 2.52 in Donea & Huerta (2003):
      
      .. math:: \int_\Omega \frac{1}{2}\frac{\bar{\kappa}}{\|\mathbf{u}\|_2^2}\left(\mathbf{u}\cdot\nabla\mathbf{w}\right)\left(\mathbf{u}\cdot\nabla\mathbf{u}\right)
      
      where :math:`\kappa` is the (isotropic) diffusivity, :math:`\mathbf{u}` is the velocity field, and :math:`\mathbf{w}` is the test function.
      
      :param ufl.TestFunction w: The test function.
      :param ufl.Function u: The velocity field
      :param ufl.Function magnitude: The magnitude of the velocity field.
      :param ufl.Function grid_pe: The grid Peclet number field.
      :returns: A UFL Form representing the streamline-upwind stabilisation term.
      :rtype: ufl.Form      
      """
      dimension = len(u)

      #FIXME: Allow the user to provide the scaling factor via the options tree.
      scaling_factor = 0.5
      k_bar = self.k_bar(magnitude, grid_pe)

      F = scaling_factor*(k_bar/(magnitude**2))*inner(dot(grad(w), u), dot(grad(u), u))*dx

      return F

   def k_bar(self, magnitude, grid_pe):
      r""" Construct :math:`\bar{\kappa}`, defined by
      
      .. math:: \bar{\kappa} = \left( \frac{1}{\tanh(\mathrm{Pe})} - \frac{1}{\mathrm{Pe}} \right) \Delta_e \|\mathbf{u}\|_2
      
      where :math:`\mathbf{u}` is the velocity field, :math:`\mathrm{Pe}` is the grid Peclet number
      
      .. math:: \mathrm{Pe} = \frac{\|\mathbf{u}\|_2\Delta_e}{2\kappa}
      
      and :math:`\Delta_e` is a measure of the cell size.
      
      :param ufl.Function magnitude: The Euclidean norm of the velocity field
      :param ufl.Function grid_pe: The grid Peclet number field.
      :returns: The term :math:`\bar{\kappa}`.
      """
      return ( (1.0/tanh(grid_pe)) - (1.0/grid_pe) ) * self.cellsize * magnitude
      
