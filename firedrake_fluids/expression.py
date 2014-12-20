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
from firedrake_fluids import LOG
import libspud

class ExpressionFromOptions:

   def __init__(self, path, t=None):

      if(libspud.have_option(path + "/constant")):
         self.val = libspud.get_option(path + "/constant")
         self.constant = True
            
      elif(libspud.have_option(path + "/python")):
         v = libspud.get_option(path + "/python")   
         self.constant = False   
         exec v
         self.val = val
         self.t = t
         
      return 
      
   def get_expression(self):
      if(self.constant):
         return Expression(self.val)
      else:
         val = self.val
         t = self.t
         # Determine the value shape by plugging in some dummy coordinate and time.
         s = val(x = [0,0,0], t=t)
         
         class PythonExpression(Expression):
            def eval(self, value, x, t=None):
               value[:] = val(x, t)
               
            if(not isinstance(s, float) and not isinstance(s, int)):
               def value_shape(self):
                  return (len(s),)

         e = PythonExpression(t=t)
         return e
