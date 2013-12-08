from firedrake import *
from random import random

def element_volume(v):
   return v.cell().volume
   
def strain_rate_tensor(v):
   return 0.5*(grad(v) + grad(v).T)

def eddy_viscosity(mesh):
   W = FunctionSpace(mesh, "CG", 1)
   w = TestFunction(W)
   eddy_viscosity = TrialFunction(W)

   if(mesh.geometry().dim() == 2):
      filter_width = element_volume(eddy_viscosity)**(1.0/2.0)
   elif(mesh.geometry().dim() == 3):
      filter_width = element_volume(eddy_viscosity)**(1.0/3.0)
   else:
      print "Dimension == 1"
      
   smagorinsky_constant = Constant(1.0)
   density = Constant(1.0)

   U = VectorFunctionSpace(mesh, "CG", 1)
   #u = Function(U)
   #for i in range(0, len(u.vector())):
   #   u.vector()[i] = random()
   u = project(Expression(('sin(x[0])', 'cos(x[0])')), U)

   #print strain_rate_tensor
   #strain_rate_tensor = inner(u, nabla_grad(u))
   shape = strain_rate_tensor(u).shape()
   second_invariant = 0.0
   for i in range(0, shape[0]):
      for j in range(0, shape[1]):
         second_invariant += 2.0*(strain_rate_tensor(u)[i,j]**2)
   second_invariant = sqrt(second_invariant)
   rhs = density*(smagorinsky_constant*filter_width)**2*second_invariant

   a = inner(w, eddy_viscosity)*dx
   L = w*rhs*dx

   solution = Function(W)
   solve(a == L, solution)

   return solution

#mesh = Mesh("mesh.xml.gz")
errors = []
for n in [2, 4, 8]:
   mesh = UnitSquare(n, n)
   smagorinsky_coefficient = 1.0*(1.0/n) # Remember to halve the Smagorinsky coefficient each time to compensate for the smaller filter width
   exact_solution = Expression('coefficient * sqrt(2.0*cos(x[0])*cos(x[0]) + 0.5*sin(x[0])*sin(x[0]) + 0.5*sin(x[0])*sin(x[0]))', coefficient = smagorinsky_coefficient)
   errors.append(errornorm(eddy_viscosity(mesh), exact_solution))

#plot(exact_solution)
#plot(solution)
#interactive()
import pylab
pylab.loglog([1.0/i for i in [2,4,8]], errors, label="Numerical error")
pylab.loglog([1.0/i for i in [2,4,8]], [0.5*(0.5**i) for i in range(0, 3)], label="First-order convergence")
pylab.xlabel(r"Characteristic element length")
pylab.ylabel(r"Error in L^2 norm")
pylab.legend()
pylab.savefig("les_convergence.png")
