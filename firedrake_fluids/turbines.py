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

from firedrake import *
      
class TopHatTurbine(Expression):

   def eval(self, value, X, K=None, coords=None, r=None):
      px = sqrt((X[0]-coords[0])**2)
      py = sqrt((X[1]-coords[1])**2)
      
      if(px <= r[0] and py <= r[1]):
         value[0] = K
      else:
         value[0] = 0
      
class BumpTurbine(Expression):

   def eval(self, value, X, K=None, coords=None, r=None):
      value[0] = K
      for i in range(2):
         p = sqrt(((X[i]-coords[i])/r[i])**2)
         if(p < 1.0):
            value[0] *= exp( 1.0 - 1.0/(1.0 - p**2) )
         else:
            value[0] *= 0

