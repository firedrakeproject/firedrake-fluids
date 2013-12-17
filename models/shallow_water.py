import sys
import libspud
import numpy

from firedrake import *

import stabilisation
import fields_calculations
import diagnostics

class VectorExpressionFromOptions(Expression):
   def __init__(self, path, t):
      self.path = path
      if(libspud.have_option(path + "/constant")):
         self.source_type = "constant"
         self.source_value = libspud.get_option(path + "/constant")
         Expression.__init__(self, code=self.source_value)
      elif(libspud.have_option(path + "/python")):
         self.source_type = "python"
         self.source_value = libspud.get_option(path + "/python")      
      elif(libspud.have_option(path + "/cpp")):
         self.source_type = "cpp"
         self.source_value = libspud.get_option(path + "/cpp")
         exec self.source_value
         Expression.__init__(self, code=val(t))
      else:
         print "No value specified"
         sys.exit(1)
   #def eval(self, values, x):
   #   if(self.constant):
   #      values[:] = self.source_value
   #   else:
   #      values[:] = self.get_values_from_python(x)
   #def get_values_from_python(self, x):
   #   s = self.source_value
   #   exec s
   #   return val(x)
   #def value_shape(self):
   #   return (2,)
      
class ScalarExpressionFromOptions(Expression):
   def __init__(self, path, t):  
      self.path = path
      if(libspud.have_option(path + "/constant")):
         self.constant = True
         self.source_value = libspud.get_option(path + "/constant")
         Expression.__init__(self, code=self.source_value)
      elif(libspud.have_option(path + "/python")):
         self.constant = False
         self.source_value = libspud.get_option(path + "/python")
      elif(libspud.have_option(path + "/cpp")):
         self.source_type = "cpp"
         self.source_value = libspud.get_option(path + "/cpp")
         exec self.source_value
         Expression.__init__(self, code=val(t))
      else:
         print "No value specified"
         sys.exit(1)
   #def eval(self, values, x):
   #   if(self.constant):
   #      values[0] = self.source_value
   #   else:
   #      values[0] = self.get_values_from_python(x)
   #def get_values_from_python(self, x):
   #   s = self.source_value
   #   exec s
   #   return val(x)

def solve_shallow_water(mesh, function_spaces):
   
   # Define function spaces
   U = function_spaces["VelocityFunctionSpace"]
   H = function_spaces["FreeSurfaceFunctionSpace"]
   W = MixedFunctionSpace([U for dim in range(dimension)] + [H])
  
   # The solution field defined on the mixed function space
   solution = Function(W)
   output_function = Function(W.sub(dimension))
  
   # Get time-stepping parameters
   T = libspud.get_option("/timestepping/finish_time")
   t = libspud.get_option("/timestepping/current_time")
   dt = libspud.get_option("/timestepping/timestep")
   if(libspud.have_option("/timestepping/steady_state")):
      steady_state_tolerance = libspud.get_option("/timestepping/steady_state/tolerance")
   else:
      steady_state_tolerance = -1000.0
       
   # Define the compulsory shallow water fields
   u_old = [Function(W.sub(dim)) for dim in range(dimension)]
   h_old = Function(W.sub(dimension))

   # Set up initial conditions
   h_initial = ScalarExpressionFromOptions(path = "/material_phase[0]/scalar_field::FreeSurfacePerturbationHeight/prognostic/initial_condition", t=t)
   h_old.interpolate(h_initial)

   expr = VectorExpressionFromOptions(path = "/material_phase[0]/vector_field::Velocity/prognostic/initial_condition", t=t)
   u_initial = [expr.code[dim] for dim in range(dimension)]   
   for dim in range(dimension):
      u_old[dim].interpolate(Expression(u_initial[dim]))
   
   # Write initial condition to file
   output_function.assign(h_old)
   output_file << output_function

   # Define trial and test functions
   functions = TrialFunctions(W)
   h = functions[-1]; u = functions[:-1]
   functions = TestFunctions(W)
   v = functions[-1]; w = functions[:-1]

   # Physical parameters
   g = libspud.get_option("/physical_parameters/gravity/magnitude")
   
   # Enable/disable terms in the shallow water equations
   if(libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/mass_term/exclude_mass_term")):
      have_momentum_mass = False
   else:
      have_momentum_mass = True
      
   if(libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/advection_term/exclude_advection_term")):
      have_momentum_advection = False
   else:
      have_momentum_advection = True
      
   if(libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/viscosity")):
      have_momentum_stress = True
      nu = libspud.get_option("/material_phase[0]/vector_field::Velocity/prognostic/viscosity")
   else:
      have_momentum_stress = False

   # Normal vector to each element facet
   n = FacetNormal(mesh.ufl_cell())

   # Mean free surface height
   h_mean = Function(W.sub(dimension))
   h_mean.interpolate(ScalarExpressionFromOptions(path = "/material_phase[0]/scalar_field::FreeSurfaceMeanHeight/prescribed/value", t=t))

   # Get any source terms for the momentum and continuity equations
   if(libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Source")):
      have_momentum_source = True
      expr = VectorExpressionFromOptions(path = "/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Source/prescribed/value", t=t)
      momentum_source = [Function(W.sub(dim)).interpolate(Expression(expr.code[dim])) for dim in range(dimension)]
   else:
      have_momentum_source = False

   if(libspud.have_option("/material_phase[0]/scalar_field::FreeSurfacePerturbationHeight/prognostic/scalar_field::Source")):
      have_continuity_source = True
      expr = ScalarExpressionFromOptions(path = "/material_phase[0]/scalar_field::FreeSurfacePerturbationHeight/prognostic/scalar_field::Source/prescribed/value", t=t)
      continuity_source = Function(W.sub(dimension)).interpolate(Expression(expr.code[0]))
   else:
      have_continuity_source = False

   # Check for any SU stabilisation
   have_momentum_su_stabilisation = libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/streamline_upwind_stabilisation")

   # Picard iteration loop settings
   max_nlit = libspud.get_option("/timestepping/nonlinear_iterations") # Maximum number of non-linear iterations allowed
   nlit_tolerance = libspud.get_option("/timestepping/nonlinear_iterations/tolerance") # Tolerance for the non-linear iteration loop
   # Non-linear fields
   u_k = [Function(W.sub(dim)) for dim in range(dimension)]
   h_k = Function(W.sub(dimension))
         
   # The time-stepping loop
   while t < T:
      t += dt
      print "\nt = %g" % t

      for dim in range(0, dimension):
         u_k[dim].assign(u_old[dim])
      h_k.assign(h_old)
         
      # Picard non-linear iteration loop
      eps = 1.0           # error measure ||u_tent - u_k||
      iter = 0            # iteration counter
      while eps > nlit_tolerance and iter < max_nlit:
         iter += 1
         
         # The collection of all the individual terms in their weak form.
         F = 0

         # Mass term
         if(have_momentum_mass):
            M_momentum = 0
            for dim in range(0, dimension):
               M_momentum += (1.0/dt)*(inner(w[dim], u[dim]) - inner(w[dim], u_old[dim]))*dx
            F += M_momentum
         
         # Advection term
         if(have_momentum_advection):
            integrate_advection_by_parts = libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/advection_term/integrate_advection_term_by_parts")
            A_momentum = 0
            if(integrate_advection_by_parts):
               raise NotImplementedError("Integration of the advection term by parts is not yet supported.")
               for dim in range(0, dimension):
                  #outflow = (dot(u_k[dim], n[dim]) + abs(dot(u_k[dim], n[dim])))/2.0
                  #expression = VectorExpressionFromOptions(path = "/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions[0]/type::dirichlet")
                  A_momentum_int = -inner(dot(u[dim], grad(w[dim])[dim]), u_k[dim])*dx - inner(w[dim]*div(u[dim]), u_k[dim])*dx
                  A_momentum_bdy = inner(w[dim], dot(u_k[dim], n[dim])*u[dim])*ds(2)
                  #A_momentum_facet = dot( outflow('+')*u[dim]('+') - outflow('-')*u[dim]('-'), jump(w[dim]))*dS
                  A_momentum += A_momentum_int + A_momentum_bdy #+ A_momentum_facet

            else:
               for dim_i in range(dimension):
                  for dim_j in range(dimension):
                     A_momentum += inner(dot(grad(u_k[dim_i])[dim_j], u[dim_j]), w[dim_i])*dx
            F += A_momentum
            
         # Viscous stress term. Note that 'nu' is the kinematic (not dynamic) viscosity.
         if(have_momentum_stress):
            K_momentum = 0
            for dim in range(dimension):
                  K_momentum += -nu*inner(grad(w[dim]), grad(u[dim]))*dx
            F -= K_momentum # Negative sign here because we are bringing the stress term over from the RHS.

         # The gradient of the height of the free surface, h
         C_momentum = 0
         for dim in range(dimension):
            C_momentum += -g*inner(w[dim], grad(h)[dim])*dx
         F -= C_momentum

         # Quadratic drag term in the momentum equation
         have_bottom_drag = libspud.have_option("/material_phase[0]/scalar_field::DragCoefficient")
         if(have_bottom_drag):
            C_D = libspud.get_option("/material_phase[0]/scalar_field::DragCoefficient/prescribed/value/constant")
            D_momentum = 0
            magnitude = 0
            for dim in range(dimension):
               magnitude += dot(u_k[dim], u_k[dim])
            magnitude = sqrt(magnitude)
            for dim in range(dimension):
               D_momentum += -inner(w[dim], (C_D*magnitude/h_mean)*u[dim])*dx
            F -= D_momentum

         # The mass term in the shallow water continuity equation 
         # (i.e. an advection equation for the free surface height, h)
         M_continuity = (1.0/dt)*(inner(v, h) - inner(v, h_old))*dx
         F += M_continuity

         # Divergence term in the shallow water continuity equation
         integrate_continuity_by_parts = libspud.have_option("/material_phase[0]/integrate_continuity_equation_by_parts")
         if(integrate_continuity_by_parts):
            Ct_continuity = 0
            # Add in any weak boundary conditions
            #expression = VectorExpressionFromOptions(path = "/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions[0]/type::dirichlet")
            for dim in range(dimension):
               Ct_continuity += - h_mean*inner(u[dim], grad(v)[dim])*dx
                                #+ h_mean * inner(u[dim], n[dim]) * v * ds
                                #+ h_mean('+')*inner(jump(v, n), avg(u[dim]))*dS
         else:
            divergence = 0
            for dim in range(dimension):
               divergence += grad(u[dim])[dim]
            Ct_continuity = h_mean*(inner(v, divergence))*dx
         F += Ct_continuity
            
         # Add in any source terms
         if(have_momentum_source):
            print "Adding momentum source..."
            for dim in range(dimension):
               F -= inner(w[dim], momentum_source[dim])*dx

         if(have_continuity_source):
            print "Adding continuity source..."
            F -= inner(v, continuity_source)*dx
            
         # Add in any SU stabilisation
         if(have_momentum_su_stabilisation):
            print "Adding momentum SU stabilisation..."
            F += stabilisation.streamline_upwind(mesh, dimension, w, u, u_k)

         # Split up the form into the LHS and RHS (i.e. the bilinear and linear forms)
         a = lhs(F)
         L = rhs(F)

         # Get all the Dirichlet boundary conditions for the Velocity field
         bcs = []
         for i in range(0, libspud.option_count("/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions/type::dirichlet")):
            if(not(libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions[%d]/type::dirichlet/apply_weakly" % i))):
               # If it's not a weak BC, then it must be a strong one.
               print "Applying Velocity BC #%d" % i
               expr = VectorExpressionFromOptions(path = ("/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions[%d]/type::dirichlet" % i), t=t)
               # Surface IDs on the domain boundary
               surface_ids = libspud.get_option("/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions[%d]/surface_ids" % i)
               for dim in range(dimension):
                  bc = DirichletBC(W.sub(dim), Expression(expr.code[dim]), surface_ids)
                  bcs.append(bc)

         # Get all the Dirichlet boundary conditions for the FreeSurfacePerturbation field
         for i in range(0, libspud.option_count("/material_phase[0]/scalar_field::FreeSurfacePerturbationHeight/prognostic/boundary_conditions/type::dirichlet")):
            if(not(libspud.have_option("/material_phase[0]/scalar_field::FreeSurfacePerturbationHeight/prognostic/boundary_conditions[%d]/type::dirichlet/apply_weakly" % i))):
               # If it's not a weak BC, then it must be a strong one.
               print "Applying FreeSurfacePerturbationHeight BC #%d" % i
               expr = ScalarExpressionFromOptions(path = ("/material_phase[0]/scalar_field::FreeSurfacePerturbationHeight/prognostic/boundary_conditions[%d]/type::dirichlet" % i), t=t)
               # Surface IDs on the domain boundary
               surface_ids = libspud.get_option("/material_phase[0]/scalar_field::FreeSurfacePerturbationHeight/prognostic/boundary_conditions[%d]/surface_ids" % i)
               bc = DirichletBC(W.sub(dimension), Expression(expr.code[0]), surface_ids)
               bcs.append(bc)
            
         # Solve the system of equations!
         solve(a == L, solution, bcs)
         fields = solution.split()
         u_tent = fields[0:-1]; h_tent = fields[-1]

         diff = abs(u_tent[0].vector().array() - u_k[0].vector().array())
         eps = numpy.linalg.norm(diff, ord=numpy.Inf)
         print 'Iteration %d: Norm=%g\n\n' % (iter, eps)

         for dim in range(0, dimension):
            u_k[dim].assign(u_tent[dim])
         h_k.assign(h_tent)
         
      # Write the solution to a file.
      print "Writing data to file..."
      output_function.assign(h_k)
      output_file << output_function
      
      # Check whether a steady-state has been reached.
      if(max(abs(h_k.vector().array() - h_old.vector().array())) <= steady_state_tolerance):
         print "Steady-state attained. Exiting the time-stepping loop..."
         break

      # Move to next time step
      for dim in range(dimension):
         u_old[dim] = project(u_k[dim], W.sub(dim))
      h_old = project(h_k, W.sub(dimension))
      print "Moving to next time level..."      

   print "Out of the time-stepping loop."
   
   # Compute the error for MMS tests
   mms = False
   if(mms):
      exact = Function(H).interpolate(Expression("sin(x[0])*sin(x[1])"))   
      test = TestFunction(H)
      trial = TrialFunction(H)
      error = Function(H)
      a = inner(test, trial)*dx
      L = inner(test, exact - project(h_old, H))*dx
      solve(a == L, error)
      print max(abs(error.vector().array()))
   

if(__name__ == "__main__"):
   # Construct the full options tree
   try:
      libspud.load_options(str(sys.argv[1]))
   except IndexError:
      print "Please provide the path to the simulation configuration file."
      sys.exit(1)

   simulation_name = libspud.get_option("/simulation_name")
   print "Initialising simulation '%s'" % simulation_name
   
   # Read in the input mesh (or construct a unit mesh)
   dimension = libspud.get_option("/geometry/dimension")
   if(libspud.have_option("/geometry/mesh/unit_mesh")):
      number_of_nodes = libspud.get_option("/geometry/mesh/unit_mesh/number_of_nodes")
      if(dimension == 1):
         mesh = UnitInterval(number_of_nodes[0])
      elif(dimension == 2):
         mesh = UnitSquareMesh(number_of_nodes[0], number_of_nodes[1])
      elif(dimension == 3):
         mesh = UnitCube(number_of_nodes[0], number_of_nodes[1], number_of_nodes[2])
      else:
         print "Unsupported dimension."
         sys.exit(1)
   elif(libspud.have_option("/geometry/mesh/from_file")):
      mesh_path = libspud.get_option("/geometry/mesh/from_file/path")
      mesh = Mesh(mesh_path)
   else:
      print "Unsupported input mesh type."
      sys.exit(1)

   # Create a dictionary containing all the function spaces
   function_spaces = {}
   for i in range(0, libspud.option_count("/function_spaces/function_space")):
      function_space_name = libspud.get_option("/function_spaces/function_space["+str(i)+"]/name")
      function_space_type = libspud.get_option("/function_spaces/function_space["+str(i)+"]/type")
      function_space_polynomial_degree = libspud.get_option("/function_spaces/function_space["+str(i)+"]/polynomial_degree")
      function_space_element_type = libspud.get_option("/function_spaces/function_space["+str(i)+"]/element_type")
      function_space_continuity = libspud.get_option("/function_spaces/function_space["+str(i)+"]/continuity")
      print "Setting up a new %s function space called %s" % (function_space_type, function_space_name)
      if(function_space_element_type == "lagrangian" and function_space_continuity == "continuous"):
         function_spaces[function_space_name] = FunctionSpace(mesh, "CG", function_space_polynomial_degree)
      else:
         print "Unknown element type and/or continuity."
         sys.exit(1)
             
   # Set up the output stream
   output_file = File("%s.pvd" % simulation_name)

   # Solve the shallow water equations!
   solve_shallow_water(mesh, function_spaces)

