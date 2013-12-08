from firedrake import *

def advect(mesh, H, u, h_old):

   w = TestFunction(H)
   h = TrialFunction(H)
  
   # Mass term
   mass_term = (1.0/Constant(dt))*inner(w, h) * dx
   
   # Advection term
   advection_term = inner(w, dot(u, grad(h))) * dx + inner(w, dot(h, div(u))) * dx
   
   # Diffusion term
   diffusion_term = Constant(0.001/1000.0)*dot(grad(w), grad(h))*dx

   a = mass_term + advection_term# - diffusion_term
   L = (1.0/Constant(dt))*inner(w, h_old) * dx

   solution = Function(H)
   solve(a == L, solution)
   return solution

