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
      
      try:
         if(libspud.have_option(path + "/constant")):
            self.val = libspud.get_option(path + "/constant")
            self.type = "constant"
               
         elif(libspud.have_option(path + "/python")):
            v = libspud.get_option(path + "/python")   
            self.type = "python"
            exec v # Make the 'val' function that the user has defined available for calling.
            self.val = val
            self.t = t
            
         elif(libspud.have_option(path + "/cpp")):
            # For C++ expressions.
            self.type = "cpp"
            v = libspud.get_option(path + "/cpp")
            exec v
            self.val = val
            self.t = t
         else:
            raise ValueError("Unknown expression type.")
      except ValueError as e:
         LOG.exception(e)
         sys.exit()
         
      return 
      
   def get_expression(self):
      """ Return a UFL Expression object, whose value is obtained from self.val. """
      try:
         if(self.type == "constant"):
            return Expression(self.val)
         elif(self.type == "cpp"):
            return Expression(code=self.val(), t=self.t)
         elif(self.type == "python"):
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
         else:
            raise ValueError("Unknown expression type: %s." % self.type)
            
      except ValueError as e:
         LOG.exception(e)
         sys.exit()
         
