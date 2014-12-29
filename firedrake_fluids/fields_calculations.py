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

def magnitude_vector(u, function_space):
   r""" Calculate the Euclidean norm of a given vector-valued Function :math:`\mathbf{u}`
   
        .. math:: \|\mathbf{u}\|_2 = \sqrt{\mathbf{u}\cdot\mathbf{u}}
        
   :param ufl.Function u: The vector field.
   :param ufl.FunctionSpace function_space: The (scalar) function space used to represent the 'magnitude' field.
   :returns: A UFL Function representing the magnitude of the given vector field.
   :rtype: ufl.Function
   """

   w = TestFunction(function_space)
   magnitude = TrialFunction(function_space)
   solution = Function(function_space)

   a = w*magnitude*dx
   L = w*sqrt(dot(u, u))*dx
   solve(a == L, solution, bcs=[])

   return solution

def grid_peclet_number(diffusivity, magnitude, function_space, cellsize):
   r""" Calculate the grid Peclet number field, given by
   
       .. math:: \mathrm{Pe} = \frac{\|\mathbf{u}\|_2\Delta x}{2\kappa}
       
       where :math:`\kappa` is the (isotropic) diffusivity, :math:`\Delta x` is the element size, and :math:`\mathbf{u}` is the velocity field.
       
   :param diffusivity: The isotropic diffusivity, which can be a constant value or a ufl.Function.
   :param ufl.Function magnitude: The magnitude of the velocity field.
   :param ufl.FunctionSpace function_space: The (scalar) function space used to represent the grid Peclet number field.
   :param ufl.CellSize cellsize: A measure of the size of the cells in the mesh.
   :returns: A UFL Function representing the grid Peclet number.
   :rtype: ufl.Function
   """

   w = TestFunction(function_space)
   grid_pe = TrialFunction(function_space)
   solution = Function(function_space)

   a = w*grid_pe*dx
   L = w*(magnitude*cellsize)/(2.0*diffusivity)*dx
   solve(a == L, solution, bcs=[])

   return solution

def steady_state(f, f_old, tolerance=1e-7):
   """ Determine whether the field 'f' (which is a Function) has reached steady-state, given a steady-state tolerance.
   
   Note that this uses the infinity norm.
   
   :param ufl.Function f: The field at time n.
   :param ufl.Function f_old: The field at time n-1.
   :param float tolerance: The steady-state tolerance.
   :returns: True if the field has reached steady-state, and False otherwise.
   :rtype: bool
    """
   
   # Take the maximum difference across all processes.
   maximum_difference = max(abs(f.vector().gather() - f_old.vector().gather()))
   if(maximum_difference <= tolerance):
      return True
   else:
      return False

