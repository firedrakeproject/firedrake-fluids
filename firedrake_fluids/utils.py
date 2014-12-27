#    Copyright (C) 2014 Imperial College London.

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

from mpi4py import MPI

def rank0():
   """ Determine whether or not the calling process has MPI rank 0.
   
   :returns: True if the calling process has MPI rank 0, otherwise returns False.
   :rtype: bool
   """
   if(MPI.COMM_WORLD.Get_rank() == 0):
      return True
   else:
      return False
