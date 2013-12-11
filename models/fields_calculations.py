from firedrake import *

def magnitude_vector(mesh, u):

   #magnitude = Function(S)
   #u_nodes = u.vector()
   #print type(u_nodes)

   #for i in range(0, len(u_nodes)):
   #   magnitude.vector()[i] = norm(u_nodes[i], norm_type='L2')
   #error = s*sqrt(dot(u, u))*dx
   #mass = inner(s, s)*dx
   #solve(mass == error, magnitude)

   FS = FunctionSpace(mesh, "CG", 1)
   w = TestFunction(FS)
   magnitude = TrialFunction(FS)
   solution = Function(FS)
   
   a = w*magnitude*dx
   s = sqrt(dot(u[0], u[0]) + dot(u[1], u[1]))
   L = w*s*dx
   solve(a == L, solution, [])
   
   return solution


