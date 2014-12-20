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

def magnitude_vector(mesh, u, function_space):
   r""" Calculate the Euclidean norm of a given vector-valued Function :math:`\mathbf{u}`. """

   w = TestFunction(function_space)
   magnitude = TrialFunction(function_space)
   solution = Function(function_space)

   a = w*magnitude*dx
   L = w*sqrt(dot(u, u))*dx
   solve(a == L, solution, bcs=[])

   return solution

def grid_peclet_number(mesh, diffusivity, magnitude, function_space, cellsize):
   r""" Calculate and return the grid Peclet number field, given by
   
       .. math:: \mathrm{Pe} = \frac{\|\mathbf{u}\|\Delta x}{2\kappa}
       
       where :math:`\kappa` is the (isotropic) diffusivity, :math:`\Delta x` is the element size, and :math:`\mathbf{u}` is the velocity field.
   """

   w = TestFunction(function_space)
   grid_pe = TrialFunction(function_space)
   solution = Function(function_space)

   a = w*grid_pe*dx
   L = w*(magnitude*cellsize)/(2.0*diffusivity)*dx
   solve(a == L, solution, bcs=[])

   return solution

def steady_state(f, f_old, tolerance=1e-7):
   """ Return True if the field 'f' (which is a Function) has reached steady-state, given a steady-state tolerance. Return False otherwise. """
   
   # Take the maximum difference across all processes.
   maximum_difference = max(abs(f.vector().gather() - f_old.vector().gather()))
   if(maximum_difference <= tolerance):
      return True
   else:
      return False

