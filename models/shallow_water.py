# Copyright 2013 Imperial College London. All rights reserved.

import sys
import libspud
import numpy
import time

backend = "firedrake"
if(backend == "firedrake"):
   from firedrake import *
else:
   from dolfin import *
   
#from firedrake_adjoint import *

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

class ShallowWater:
   
   def __init__(self, path):
      """ Initialise a new shallow water simulation using an options file stored at the location given by 'path'. """
   
      print "Initialising simulation..."
      
      # Remove any stored options.
      libspud.clear_options()

      # Load the options from the options tree.
      libspud.load_options(path)
     
      # Populate the options dictionary
      self.populate_options()
        
      # Read in the input mesh (or construct a unit mesh)
      dimension = self.options["dimension"]
      if(libspud.have_option("/geometry/mesh/unit_mesh")):
         number_of_nodes = libspud.get_option("/geometry/mesh/unit_mesh/number_of_nodes")
         if(dimension == 1):
            self.mesh = UnitIntervalMesh(number_of_nodes[0])
         elif(dimension == 2):
            self.mesh = UnitSquareMesh(number_of_nodes[0], number_of_nodes[1])
         elif(dimension == 3):
            self.mesh = UnitCubeMesh(number_of_nodes[0], number_of_nodes[1], number_of_nodes[2])
         else:
            print "Unsupported dimension."
            sys.exit(1)
      elif(libspud.have_option("/geometry/mesh/from_file")):
         path_to_config = os.path.dirname(os.path.abspath(path))
         path_to_mesh = libspud.get_option("/geometry/mesh/from_file/path") # This is the path relative to the directory where the configuration file is stored.
         self.mesh = Mesh(os.path.join(path_to_config, path_to_mesh))
      else:
         print "Unsupported input mesh type."
         sys.exit(1)

      # Create a dictionary containing all the function spaces
      self.function_spaces = {}
      for i in range(0, libspud.option_count("/function_spaces/function_space")):
         function_space_name = libspud.get_option("/function_spaces/function_space["+str(i)+"]/name")
         function_space_type = libspud.get_option("/function_spaces/function_space["+str(i)+"]/type")
         function_space_polynomial_degree = libspud.get_option("/function_spaces/function_space["+str(i)+"]/polynomial_degree")
         function_space_element_type = libspud.get_option("/function_spaces/function_space["+str(i)+"]/element_type")
         function_space_continuity = libspud.get_option("/function_spaces/function_space["+str(i)+"]/continuity")
         print "Setting up a new %s function space called %s" % (function_space_type, function_space_name)
         if(function_space_element_type == "lagrangian" and function_space_continuity == "continuous"):
            self.function_spaces[function_space_name] = FunctionSpace(self.mesh, "CG", function_space_polynomial_degree)
         else:
            print "Unknown element type and/or continuity."
            sys.exit(1)
                
      # Define the mixed function space
      U = self.function_spaces["VelocityFunctionSpace"]
      H = self.function_spaces["FreeSurfaceFunctionSpace"]
      self.W = MixedFunctionSpace([U for dim in range(dimension)] + [H])
     
      # Define trial and test functions
      functions = TrialFunctions(self.W)
      self.h = functions[-1]; self.u = functions[:-1]
      functions = TestFunctions(self.W)
      self.v = functions[-1]; self.w = functions[:-1]

      # Normal vector to each element facet
      if(backend == "firedrake"):
         self.n = FacetNormal(self.mesh.ufl_cell())
      else:
         self.n = FacetNormal(self.mesh)
      
      # The solution field defined on the mixed function space
      self.solution = Function(self.W)
      self.output_function = Function(self.W.sub(dimension), name="FreeSurfacePerturbation")
     
      # Define the compulsory shallow water fields
      self.u_old = [Function(self.W.sub(dim)) for dim in range(dimension)]
      self.h_old = Function(self.W.sub(dimension))

      # Non-linear fields
      self.u_k = [Function(self.W.sub(dim)) for dim in range(dimension)]
      self.h_k = Function(self.W.sub(dimension))
      
      # Set up initial conditions
      h_initial = ScalarExpressionFromOptions(path = "/material_phase[0]/scalar_field::FreeSurfacePerturbation/prognostic/initial_condition", t=self.options["t"])
      self.h_old.interpolate(h_initial)

      expr = VectorExpressionFromOptions(path = "/material_phase[0]/vector_field::Velocity/prognostic/initial_condition", t=self.options["t"])
      u_initial = [expr.code[dim] for dim in range(dimension)]   
      for dim in range(dimension):
         self.u_old[dim].interpolate(Expression(u_initial[dim]))
   
      # Mean free surface height
      self.h_mean = Function(self.W.sub(dimension))
      self.h_mean.interpolate(ScalarExpressionFromOptions(path = "/material_phase[0]/scalar_field::FreeSurfaceMean/prescribed/value", t=self.options["t"]))

      # Set up the output stream
      self.output_file = File("%s.pvd" % self.options["simulation_name"])
      
      # Write initial condition to file
      self.output_function.assign(self.h_old)
      self.output_file << self.output_function
            
      return
      
   def populate_options(self):
      """ Add simulation options related to the shallow water model to a dictionary object. """
      # A dictionary storing all the options
      self.options = {}
      
      self.options["simulation_name"] = libspud.get_option("/simulation_name")
      self.options["dimension"] = libspud.get_option("/geometry/dimension")
   
      # Time-stepping parameters
      self.options["T"] = libspud.get_option("/timestepping/finish_time")
      self.options["t"] = libspud.get_option("/timestepping/current_time")
      self.options["dt"] = libspud.get_option("/timestepping/timestep")
      if(libspud.have_option("/timestepping/steady_state")):
         self.options["steady_state_tolerance"] = libspud.get_option("/timestepping/steady_state/tolerance")
      else:
         self.options["steady_state_tolerance"] = -1000.0

      # Picard iteration loop settings
      self.options["max_nlit"] = libspud.get_option("/timestepping/nonlinear_iterations") # Maximum number of non-linear iterations allowed
      self.options["nlit_tolerance"] = libspud.get_option("/timestepping/nonlinear_iterations/tolerance") # Tolerance for the non-linear iteration loop
       
      # Physical parameters
      self.options["g_magnitude"] = libspud.get_option("/physical_parameters/gravity/magnitude")
      
      # Enable/disable terms in the shallow water equations
      if(libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/mass_term/exclude_mass_term")):
         self.options["have_momentum_mass"] = False
      else:
         self.options["have_momentum_mass"] = True
         
      if(libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/advection_term/exclude_advection_term")):
         self.options["have_momentum_advection"] = False
      else:
         self.options["have_momentum_advection"] = True
         
      if(libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/viscosity")):
         self.options["have_momentum_stress"] = True
         self.options["nu"] = libspud.get_option("/material_phase[0]/vector_field::Velocity/prognostic/viscosity")
      else:
         self.options["have_momentum_stress"] = False
         
      # Source terms for the momentum and continuity equations
      self.options["have_momentum_source"] = libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Source")
      self.options["have_continuity_source"] = libspud.have_option("/material_phase[0]/scalar_field::FreeSurfacePerturbation/prognostic/scalar_field::Source")

      # Check for any SU stabilisation
      self.options["have_su_stabilisation"] = libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/spatial_discretisation/continuous_galerkin/streamline_upwind_stabilisation")

      self.options["have_bottom_drag"] = libspud.have_option("/material_phase[0]/scalar_field::DragCoefficient")

      self.options["integrate_continuity_by_parts"] = libspud.have_option("/material_phase[0]/integrate_continuity_equation_by_parts")

      return
      
   def run(self):
      """ Execute the time-stepping loop. """
   
      T = self.options["T"]
      t = self.options["t"]
      dt = self.options["dt"]
      dimension = self.options["dimension"]
      nlit_tolerance = self.options["nlit_tolerance"]
      max_nlit = self.options["max_nlit"]
      g_magnitude = self.options["g_magnitude"]
            
      # The time-stepping loop
      while t < T:
         t += dt
         print "\nt = %g" % t

         for dim in range(0, dimension):
            self.u_k[dim].assign(self.u_old[dim])
         self.h_k.assign(self.h_old)
            
         # Picard non-linear iteration loop
         eps = 1.0           # error measure ||u_tent - u_k||
         iter = 0            # iteration counter
         while eps > nlit_tolerance and iter < max_nlit:
            iter += 1
            
            # The collection of all the individual terms in their weak form.
            F = 0

            # Mass term
            if(self.options["have_momentum_mass"]):
               M_momentum = 0
               for dim in range(0, dimension):
                  M_momentum += (1.0/dt)*(inner(self.w[dim], self.u[dim]) - inner(self.w[dim], self.u_old[dim]))*dx
               F += M_momentum
            
            # Advection term
            if(self.options["have_momentum_advection"]):
               A_momentum = 0
               for dim_i in range(dimension):
                  for dim_j in range(dimension):
                     A_momentum += inner(dot(grad(self.u_k[dim_i])[dim_j], self.u[dim_j]), self.w[dim_i])*dx
               F += A_momentum
               
            # Viscous stress term. Note that 'nu' is the kinematic (not dynamic) viscosity.
            if(self.options["have_momentum_stress"]):
               K_momentum = 0
               for dim in range(dimension):
                     K_momentum += -self.options["nu"]*inner(grad(self.w[dim]), grad(self.u[dim]))*dx
               F -= K_momentum # Negative sign here because we are bringing the stress term over from the RHS.

            # The gradient of the height of the free surface, h
            C_momentum = 0
            for dim in range(dimension):
               C_momentum += -g_magnitude*inner(self.w[dim], grad(self.h)[dim])*dx
            F -= C_momentum

            # Quadratic drag term in the momentum equation
            
            if(self.options["have_bottom_drag"]):
               print "Adding bottom drag..."
               C_D = libspud.get_option("/material_phase[0]/scalar_field::DragCoefficient/prescribed/value/constant")
               D_momentum = 0
               magnitude = 0
               for dim in range(dimension):
                  magnitude += dot(self.u_k[dim], self.u_k[dim])
               magnitude = sqrt(magnitude)
               for dim in range(dimension):
                  D_momentum += -inner(self.w[dim], (C_D*magnitude/self.h_mean)*self.u[dim])*dx
               F -= D_momentum

            # The mass term in the shallow water continuity equation 
            # (i.e. an advection equation for the free surface height, h)
            M_continuity = (1.0/dt)*(inner(self.v, self.h) - inner(self.v, self.h_old))*dx
            F += M_continuity

            # Divergence term in the shallow water continuity equation
            if(self.options["integrate_continuity_by_parts"]):
               Ct_continuity = 0

               for dim in range(dimension):
                  Ct_continuity += - self.h_mean*inner(self.u[dim], grad(self.v)[dim])*dx
                                   #+ inner(jump(v, n)[dim], avg(u[dim]))*dS
                                 
               # Add in the surface integrals, but check to see if a no normal flow boundary condition needs to be applied weakly here.
               boundary_markers = set(self.mesh.exterior_facets.markers)
               for marker in boundary_markers:
                  marker = int(marker) # ds() will not accept markers of type 'numpy.int32', so convert it to type 'int' here.
                  no_normal_flow = False
                  for i in range(0, libspud.option_count("/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions")):
                     if(libspud.have_option("/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions[%d]/type::no_normal_flow" % i)):
                        if(marker in libspud.get_option("/material_phase[0]/vector_field::Velocity/prognostic/boundary_conditions[%d]/surface_ids" % i)):
                           no_normal_flow = True
                  if(not no_normal_flow):
                     for dim in range(dimension):
                        Ct_continuity += self.h_mean*inner(self.u[dim], self.n[dim]) * self.v * ds(int(marker))

            else:
               divergence = 0
               for dim in range(dimension):
                  divergence += grad(self.u[dim])[dim]
               Ct_continuity = self.h_mean*(inner(self.v, divergence))*dx
            F += Ct_continuity
               
            # Add in any source terms
            if(self.options["have_momentum_source"]):
               expr = VectorExpressionFromOptions(path = "/material_phase[0]/vector_field::Velocity/prognostic/vector_field::Source/prescribed/value", t=t)
               momentum_source = [Function(self.W.sub(dim)).interpolate(Expression(expr.code[dim])) for dim in range(dimension)]
               print "Adding momentum source..."
               for dim in range(dimension):
                  F -= inner(self.w[dim], momentum_source[dim])*dx

            if(self.options["have_continuity_source"]):
               expr = ScalarExpressionFromOptions(path = "/material_phase[0]/scalar_field::FreeSurfacePerturbation/prognostic/scalar_field::Source/prescribed/value", t=t)
               continuity_source = Function(self.W.sub(dimension)).interpolate(Expression(expr.code[0]))
               print "Adding continuity source..."
               F -= inner(self.v, continuity_source)*dx
               
            # Add in any SU stabilisation
            if(self.options["have_su_stabilisation"]):
               print "Adding momentum SU stabilisation..."
               F += stabilisation.streamline_upwind(self.mesh, dimension, self.w, self.u, self.u_k)

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
                     bc = DirichletBC(self.W.sub(dim), Expression(expr.code[dim]), surface_ids)
                     bcs.append(bc)

            # Get all the Dirichlet boundary conditions for the FreeSurfacePerturbation field
            for i in range(0, libspud.option_count("/material_phase[0]/scalar_field::FreeSurfacePerturbation/prognostic/boundary_conditions/type::dirichlet")):
               if(not(libspud.have_option("/material_phase[0]/scalar_field::FreeSurfacePerturbation/prognostic/boundary_conditions[%d]/type::dirichlet/apply_weakly" % i))):
                  # If it's not a weak BC, then it must be a strong one.
                  print "Applying FreeSurfacePerturbation BC #%d" % i
                  expr = ScalarExpressionFromOptions(path = ("/material_phase[0]/scalar_field::FreeSurfacePerturbation/prognostic/boundary_conditions[%d]/type::dirichlet" % i), t=t)
                  # Surface IDs on the domain boundary
                  surface_ids = libspud.get_option("/material_phase[0]/scalar_field::FreeSurfacePerturbation/prognostic/boundary_conditions[%d]/surface_ids" % i)
                  bc = DirichletBC(self.W.sub(dimension), Expression(expr.code[0]), surface_ids)
                  bcs.append(bc)
               
            # Split up the form into the LHS and RHS (i.e. the bilinear and linear forms)
            a = lhs(F)
            L = rhs(F)
            
            start = time.time()
            A = assemble(a)
            b = assemble(L)
            for bc in bcs:
               bc.apply(A)
               bc.apply(b)
            end = time.time()
            difference = end - start
            print "Tictoc 1 = %f" % difference
                      
            # Solve the system of equations!
            solution = Function(self.W)
            start = time.time()
            #solve(A, solution, b, solver_parameters={'ksp_monitor':True})
            solve(a == L, solution, bcs=bcs, solver_parameters={'ksp_monitor': True})
            end = time.time()
            difference = end - start
            print "Tictoc 2 = %f" % difference
            fields = solution.split()
            u_tent = fields[0:-1]; h_tent = fields[-1]
            
            diff = abs(u_tent[0].vector().array() - self.u_k[0].vector().array())
            eps = numpy.linalg.norm(diff, ord=numpy.Inf)
            print 'Iteration %d: Norm=%g\n\n' % (iter, eps)

            for dim in range(0, dimension):
               self.u_k[dim].assign(u_tent[dim])
            self.h_k.assign(h_tent)
            
         # Write the solution to a file.
         print "Writing data to file..."
         self.output_function.assign(self.h_k)
         self.output_file << self.output_function
         
         # Check whether a steady-state has been reached.
         if(max(abs(self.h_k.vector().array() - self.h_old.vector().array())) <= self.options["steady_state_tolerance"]):
            print "Steady-state attained. Exiting the time-stepping loop..."
            break

         # Move to next time step
         for dim in range(dimension):
            self.u_old[dim] = project(self.u_k[dim], self.W.sub(dim))
         self.h_old = project(self.h_k, self.W.sub(dimension))
         print "Moving to next time level..."      

      print "Out of the time-stepping loop."
      
      # Compute the error for MMS tests
      mms = True
      if(mms):
         H = self.function_spaces["FreeSurfaceFunctionSpace"]
         exact = Function(H).interpolate(Expression("sin(x[0])*sin(x[1])"))   
         test = TestFunction(H)
         trial = TrialFunction(H)
         error = Function(H)
         a = inner(test, trial)*dx
         L = inner(test, exact - project(self.h_old, H))*dx
         solve(a == L, error)
         f = Function(H)
         f.interpolate(Expression("sin(x[0])*sin(x[1])"))
         print sqrt(assemble(dot(self.h_old - f, self.h_old - f) * dx))
      
      return

if(__name__ == "__main__"):

   try:
      # Set up a shallow water simulation.
      sw = ShallowWater(path = sys.argv[1])
   except IndexError:
      print "Please provide the path to the simulation configuration file."
      sys.exit(1)

   # Solve the shallow water equations!
   sw.run()
   

