from firedrake import *

import fields_calculations

def streamline_upwind(mesh, v, s, u_nl):

   h = CellSize(mesh)
   scaling_factor = 0.5
   u_magnitude = fields_calculations.magnitude_vector(u_nl)

   u_nodes = u_magnitude.vector()
   for i in range(0, len(u_nodes)):
      if(u_nodes[i] < DOLFIN_EPS):
         u_nodes[i] = DOLFIN_EPS

   #grid_pe = (u_magnitude*h)/(2.0*k)

   #k_bar = ( (1.0/tanh(grid_pe)) - (1.0/grid_pe) ) * k * grid_pe

   k_bar = scaling_factor*h*u_magnitude

   form = (k_bar/(u_magnitude**2))*inner(dot(u_nl, grad(v)), dot(u_nl, grad(s)))*dx

   return form

