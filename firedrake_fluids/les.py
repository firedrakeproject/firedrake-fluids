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
from firedrake_fluids import LOG

class LES:
   """ Module containing models for Large Eddy Simulation (LES). """

   def __init__(self, mesh, function_space):
      self.mesh = mesh
      self.function_space = function_space
      return
      
   def strain_rate_tensor(self, u):
      r""" Return the UFL of the strain rate tensor :math:`\mathbb{S}`, defined by
      
      .. math:: \mathbb{S} = \frac{1}{2}\left(\nabla\mathbf{u} + \nabla\mathbf{u}^{\mathrm{T}}\right)
      
      where :math:`\mathbf{u}` is the velocity field, and the superscript T denotes the transpose.
      
      :param ufl.Function u: The velocity field.
      :returns: A UFL Form object representing the strain rate tensor.
      :rtype: ufl.Form      
      """
      S = 0.5*(grad(u) + grad(u).T)
      return S

   def eddy_viscosity(self, u, density, smagorinsky_coefficient):
      r""" Define a Form representing the eddy viscosity
      
      .. math:: \left(C_s \Delta_e\right)^2\|\mathbb{S}\|
      
      where :math:`C_s` is the Smagorinsky coefficient, :math:`\Delta_e` is a measure of the cell size, 
      and :math:`\|\mathbb{S}\|` is the modulus of the strain rate tensor :math:`\mathbb{S}`.

      :param ufl.Function u: The velocity field.
      :param density: The density field. This can be a constant value, or a ufl.Function.
      :param float smagorinsky_coefficient: The Smagorinsky coefficient. Typical values are in the range [0.1, 0.2].
      :returns: A tuple containing the LHS and RHS of the UFL Form representing the eddy viscosity.
      :rtype: tuple
      """
      dimension = len(u)
      w = TestFunction(self.function_space)
      eddy_viscosity = TrialFunction(self.function_space)

      filter_width = CellVolume(self.mesh)**(1.0/dimension)

      S = self.strain_rate_tensor(u)
      second_invariant = 0.0
      for i in range(0, dimension):
         for j in range(0, dimension):
            second_invariant += 2.0*(S[i,j]**2)
            
      second_invariant = sqrt(second_invariant)
      rhs = density*(smagorinsky_coefficient*filter_width)**2*second_invariant

      lhs = inner(w, eddy_viscosity)*dx
      rhs = inner(w, rhs)*dx
      
      return lhs, rhs

