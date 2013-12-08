from firedrake import *

def courant_number(mesh, u, dt):
   h = CellSize(mesh)
   magnitude = magnitude_vector(u)
   
   courant_number = project((magnitude*Constant(dt))/h)
   return courant_number
   
def divergence(H, vector):

   w = TestFunction(H)
   u = TrialFunction(H)

   a = inner(w, u)*dx
   L = w*div(vector)*dx

   solution = Function(H)
   solve(a == L, solution)
   return solution

def shallow_water_continuity(H, h_old, h_new, u_new):

   w = TestFunction(H)
   u = TrialFunction(H)

   a = inner(w, u)*dx
   L = w*((h_new - h_old)/Constant(dt))*dx + inner(w, dot(grad(h_new), u_new))*dx + inner(w, div(u_new)*h_new)*dx

   solution = Function(H)
   solve(a == L, solution)
   return solution

