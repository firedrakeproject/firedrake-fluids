from firedrake import *

def magnitude_vector(u):

   #magnitude = Function(S)
   #u_nodes = u.vector()
   #print type(u_nodes)

   #for i in range(0, len(u_nodes)):
   #   magnitude.vector()[i] = norm(u_nodes[i], norm_type='L2')
   #error = s*sqrt(dot(u, u))*dx
   #mass = inner(s, s)*dx
   #solve(mass == error, magnitude)

   magnitude = project(sqrt(dot(u, u)))
   return magnitude


