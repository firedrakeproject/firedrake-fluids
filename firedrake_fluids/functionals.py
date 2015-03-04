#    Copyright (C) 2015 Imperial College London.

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
from firedrake_adjoint import *

class PowerFunctional:

   def __init__(self):
      pass
      
   def Jm(self, u, drag_coefficient, density):
      return self.power(u, drag_coefficient, density)*dx
   
   def power(self, u, drag_coefficient, density):
      """ Returns the form of the power expression. """
      return density*drag_coefficient*dot(u, u)**1.5
      
