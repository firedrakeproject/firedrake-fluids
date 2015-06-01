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
from random import random
import numpy
import sys

class KEpsilon:

   def __init__(self, mesh, function_space):
      self.mesh = mesh
      self.P1 = FunctionSpace(mesh, "CG", 1)
      self.function_space = FunctionSpace(mesh, "CG", 2)

      self.w_k = TestFunction(self.function_space)
      self.k = TrialFunction(self.function_space)
      self.k_old = Function(self.function_space)
 
      self.w_eps = TestFunction(self.function_space)
      self.eps = TrialFunction(self.function_space)
      self.eps_old = Function(self.function_space)
           
      self.w_nu_t = TestFunction(self.function_space)
      self.nu_t = TrialFunction(self.function_space)
      self.nu_t_old = Function(self.function_space)
      
      self.k_old.interpolate(Expression(1.0e-9))
      self.eps_old.interpolate(Expression(1.0e-9))
      self.nu_t_old.interpolate(Expression(1.0e-9))  
      
      self.bcs = [DirichletBC(self.function_space, Expression("0.0"), (3,4,5,6))]
      
      return
      
   def strain_rate_tensor(self, u):
      S = 0.5*(grad(u) + grad(u).T)
      return S

   def eddy_viscosity(self, u, dt):

      # Constants
      C1 = 1.44
      C2 = 1.92
      C3 = -0.33
      Cnu = 0.09
      sigma_k = 1.0
      sigma_eps = 1.3

      dimension = len(u)

      # Compute modulus of strain rate tensor
      S = self.strain_rate_tensor(u)
      second_invariant = 0.0
      for i in range(0, dimension):
         for j in range(0, dimension):
            second_invariant += 2.0*(S[i,j]**2)
      second_invariant = sqrt(second_invariant)
      
      # Mass for k
      M_k = (1.0/dt)*(inner(self.w_k, self.k) - inner(self.w_k, self.k_old))*dx
      
      # Advection for k
      A_k = inner(self.w_k, div(u*self.k))*dx
      
      # Diffusion for k
      D_k = -inner(grad(self.w_k), (1e-6 + self.nu_t_old/sigma_k)*grad(self.k))*dx
      
      # Production for k
      P_k = inner(self.w_k, self.nu_t_old*inner(grad(u), grad(u)))*dx
      
      # Absorption for k
      ABS_k = -inner(self.w_k, self.eps_old)*dx

      # The full weak form of the equation for k
      F_k = M_k + A_k - D_k - P_k - ABS_k

      solution = Function(self.function_space)
      solve(lhs(F_k) == rhs(F_k), solution, bcs=self.bcs)
      self.k_old.assign(solution)

      nodes = self.k_old.vector()
      near_zero = numpy.array([1.0e-9 for i in range(len(nodes))])
      nodes.set_local(numpy.maximum(nodes.array(), near_zero))
      
      File("k.pvd") << self.k_old
      

      # Mass for epsilon
      M_eps = (1.0/dt)*(inner(self.w_eps, self.eps) - inner(self.w_eps, self.eps_old))*dx
      
      # Advection for epsilon
      A_eps = inner(self.w_eps, div(u*self.eps))*dx
      
      # Diffusion for epsilon
      D_eps = -inner(grad(self.w_eps), (1e-6 + self.nu_t_old/sigma_eps)*grad(self.eps))*dx
      
      # Production for epsilon
      P_eps = C1*(self.eps/self.k_old)*inner(self.w_eps, self.nu_t_old*inner(grad(u), grad(u)))*dx
      
      # Absorption for epsilon
      ABS_eps = -C2*(self.eps_old**2/self.k_old)*inner(self.w_eps, self.eps)*dx

      # The full weak form of the equation for epsilon
      F_eps = M_eps + A_eps - D_eps - P_eps - ABS_eps
      
      solution = Function(self.function_space)
      solve(lhs(F_eps) == rhs(F_eps), solution, bcs=self.bcs)
      self.eps_old.assign(solution)

      nodes = self.eps_old.vector()
      near_zero = numpy.array([1.0e-9 for i in range(len(nodes))])
      nodes.set_local(numpy.maximum(nodes.array(), near_zero))
      
      File("eps.pvd") << self.eps_old
      
      # Eddy viscosity
      E = inner(self.w_nu_t, self.nu_t)*dx - inner(self.w_nu_t, Cnu*((self.k_old**2)/self.eps_old))*dx
      solution = Function(self.function_space)
      solve(lhs(E) == rhs(E), solution, bcs=self.bcs)
      self.nu_t_old.assign(solution)
      
      nodes = self.nu_t_old.vector()
      near_zero = numpy.array([1.0e-9 for i in range(len(nodes))])
      nodes.set_local(numpy.maximum(nodes.array(), near_zero))
         
      File("nu_t.pvd") << self.nu_t_old
      
      import sys
      #sys.exit(1)
    
      return project(self.nu_t_old, self.P1)

