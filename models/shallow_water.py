# Copyright 2013 Imperial College London. All rights reserved.

import sys, os
import libspud
import numpy

from firedrake import *

import fields_calculations
import diagnostics
from stabilisation import *
from les import *

# FIXME: the detectors module currently relies on vtktools.
# Temporarily wrap this in a try-except block in case vtktools doesn't exist.
try:
   import detectors
except:
   print "Could not import the detectors module"

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
      elif(libspud.have_option("/geometry/mesh/interval_mesh")):
         L = libspud.get_option("/geometry/mesh/interval_mesh/length")
         n = libspud.get_option("/geometry/mesh/interval_mesh/number_of_cells")
         self.mesh = IntervalMesh(n, L)
      elif(libspud.have_option("/geometry/mesh/from_file")):
         path_to_config = os.path.dirname(os.path.abspath(path))
         # This is the path relative to the directory where the configuration file is stored.
         path_to_mesh = libspud.get_option("/geometry/mesh/from_file/relative_path") 
         self.mesh = Mesh(os.path.join(path_to_config, path_to_mesh))
      else:
         print "Unsupported input mesh type."
         sys.exit(1)

      # Create a dictionary containing all the function spaces
      self.function_spaces = {}
      for i in range(0, libspud.option_count("/function_spaces/function_space")):
         name = libspud.get_option("/function_spaces/function_space["+str(i)+"]/name")
         family = libspud.get_option("/function_spaces/function_space["+str(i)+"]/family")
         degree = libspud.get_option("/function_spaces/function_space["+str(i)+"]/degree")
         print "Setting up a new %s function space of degree %d called %s" % (family, degree, name)
         if(family == "Continuous Lagrange"):
            self.function_spaces[name] = FunctionSpace(self.mesh, "CG", degree)
         elif(family == "Discontinuous Lagrange"):
            self.function_spaces[name] = FunctionSpace(self.mesh, "DG", degree)
         else:
            print "Unknown element family: %s." % family
            sys.exit(1)
            
      # Define the coordinate function space as P1.
      self.function_spaces["CoordinateFunctionSpace"] = FunctionSpace(self.mesh, "CG", 1)
                
      # Define the mixed function space
      U = self.function_spaces["VelocityFunctionSpace"]
      H = self.function_spaces["FreeSurfaceFunctionSpace"]
      self.W = MixedFunctionSpace([U for dim in range(dimension)] + [H])
     
      # The solution field defined on the mixed function space
      self.solution = Function(self.W)
      # These are like the TrialFunctions, but are just regular Functions here because we want to solve a non-linear problem
      functions = split(self.solution)
      self.u = list(functions[:-1]); self.h = functions[-1]

      # Get the test functions
      test_functions = TestFunctions(self.W)
      self.w = test_functions[:-1]; self.v = test_functions[-1]

      # Normal vector to each element facet
      self.n = FacetNormal(self.mesh)
      
      # Set up initial conditions
      # FIXME: Subclassing of the Expression class needs to be DOLFIN compatible. The current method used here is a hack.
      h_initial = ScalarExpressionFromOptions(path = "/core_fields/scalar_field::FreeSurfacePerturbation/initial_condition", t=self.options["t"])
      expr = VectorExpressionFromOptions(path = "/core_fields/vector_field::Velocity/initial_condition", t=self.options["t"])
      u_initial = [expr.code[dim] for dim in range(dimension)]

      # Define the compulsory shallow water fields
      codes = tuple(u_initial) + (h_initial.code[0],)
      self.solution_old = Function(self.W).interpolate(Expression(codes))
      functions_old = split(self.solution_old)
      self.u_old = list(functions_old[:-1]); self.h_old = functions_old[-1]
      
      # The solution should first hold the initial condition.
      self.solution.assign(self.solution_old)
      
      # Mean free surface height
      self.h_mean = Function(self.W.sub(dimension))
      self.h_mean.interpolate(ScalarExpressionFromOptions(path = "/core_fields/scalar_field::FreeSurfaceMean/value", t=self.options["t"]))

      # Set up the functions used to write fields to file.
      self.output_function = [Function(self.W.sub(dimension), name="Velocity_%d" % dim) for dim in range(dimension)] + [Function(self.W.sub(dimension), name="FreeSurfacePerturbation")]
      
      # Set up the output stream
      self.output_file = [File("%s_Velocity_%d.pvd" % (self.options["simulation_name"], dim)) for dim in range(dimension)] + [File("%s_FreeSurfacePerturbation.pvd" % self.options["simulation_name"])]
      
      # Write initial conditions to file
      for dim in range(dimension):
         projected = project(self.solution_old.split()[dim], self.W.sub(dimension)) # Projects the Velocity to the FreeSurfaceFunctionSpace
         self.output_function[dim].assign(projected)
         self.output_file[dim] << self.output_function[dim]
      self.output_function[dimension].assign(self.solution_old.split()[dimension])
      self.output_file[dimension] << self.output_function[dimension]
    
      # Initialise detectors
      if(self.options["have_detectors"]):
         self.detectors = detectors.Detectors(locations_file_name = "detectors.xy", 
                                              values_file_name = "detectors.dat", 
                                              fields = ["FreeSurfacePerturbation", "Velocity_0", "Velocity_1"])
         self.detectors.write(self.options["simulation_name"], 0.0, self.options["dt"])
         
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
       
      # Physical parameters
      self.options["g_magnitude"] = libspud.get_option("/physical_parameters/gravity/magnitude")
      
      # Enable/disable terms in the shallow water equations
      if(libspud.have_option("/equations/momentum_equation/mass_term/exclude_mass_term")):
         self.options["have_momentum_mass"] = False
      else:
         self.options["have_momentum_mass"] = True
         
      if(libspud.have_option("/equations/momentum_equation/advection_term/exclude_advection_term")):
         self.options["have_momentum_advection"] = False
      else:
         self.options["have_momentum_advection"] = True
         
      self.options["have_momentum_stress"] = libspud.have_option("/equations/stress_term")
         
      # Source terms for the momentum and continuity equations
      self.options["have_momentum_source"] = libspud.have_option("/equations/momentum_equation/source_term")
      self.options["have_continuity_source"] = libspud.have_option("/equations/continuity_equation/source_term")

      # Check for any SU stabilisation
      self.options["have_su_stabilisation"] = libspud.have_option("/equations/momentum_equation/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_stabilisation")
      
      # Turbulence parameterisation
      self.options["have_turbulence_parameterisation"] = libspud.have_option("/equations/momentum_equation/turbulence_parameterisation")

      self.options["have_bottom_drag"] = libspud.have_option("/equations/momentum_equation/drag_term")
         
      self.options["integrate_continuity_equation_by_parts"] = libspud.have_option("/equations/continuity_equation/integrate_by_parts")

      self.options["have_detectors"] = libspud.have_option("/io/detectors/")
      
      return
      
   def run(self):
      """ Execute the time-stepping loop. """
   
      T = self.options["T"]
      t = self.options["t"]
      dt = self.options["dt"]
      dimension = self.options["dimension"]
      g_magnitude = self.options["g_magnitude"]
      H = self.h_mean + self.h # The total height of the free surface.
      
      P1 = self.function_spaces["CoordinateFunctionSpace"]
      cellsize = CellSize(self.mesh)
            
      t = dt
            
      # The time-stepping loop
      EPSILON = 1.0e-14
      while t <= T + EPSILON: # A small value EPSILON is added here in case of round-off error.
         print "\nt = %g" % t
            
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
                  A_momentum += inner(dot(grad(self.u[dim_i])[dim_j], self.u[dim_j]), self.w[dim_i])*dx
            F += A_momentum
            
         # Viscous stress term. Note that the viscosity is kinematic (not dynamic).
         if(self.options["have_momentum_stress"]):
            viscosity = ScalarExpressionFromOptions(path = "/equations/momentum_equation/stress_term/scalar_field::Viscosity/value", t=self.options["t"])
            viscosity = Function(self.W.sub(0)).interpolate(Expression(viscosity)) # Background viscosity
            if(self.options["have_turbulence_parameterisation"]):
               base_option_path = "/equations/momentum_equation/turbulence_parameterisation"
               # Add on eddy viscosity, if necessary
               if(libspud.have_option(base_option_path + "/les")):
                  les = LES(self.mesh, self.W.sub(0))
                  density = Function(self.W.sub(0)).interpolate(Expression("1.0")) # We divide through by density in the momentum equation, so just set this to 1.0 for now.
                  smagorinsky_coefficient = Function(self.W.sub(0)).interpolate(Expression(libspud.get_option(base_option_path + "/les/smagorinsky/smagorinsky_coefficient")))
                  filter_width = Function(self.W.sub(0)).interpolate(Expression(libspud.get_option(base_option_path + "/les/smagorinsky/filter_width"))) # FIXME: Remove this when CellSize is supported in Firedrake.
                  eddy_viscosity = les.eddy_viscosity(self.u, density, smagorinsky_coefficient, filter_width)
                  
               viscosity += eddy_viscosity
               
            K_momentum = 0
            
            # Compute the divergence of the velocity field.
            velocity_divergence = sum([grad(self.u[dim])[dim] for dim in range(dimension)])
            
            # Perform a double dot product of the stress tensor and grad(w).
            # tau = grad(u) + transpose(grad(u)) - (2/3)*div(u)
            for dim_i in range(dimension):
               for dim_j in range(dimension):
                  K_momentum += -viscosity*inner(grad(self.u[dim_i])[dim_j] + grad(self.u[dim_j])[dim_i], grad(self.w[dim_i])[dim_j])*dx
               K_momentum += viscosity*(2.0/3.0)*inner(velocity_divergence, grad(self.w[dim_i])[dim_i])*dx
               
            F -= K_momentum # Negative sign here because we are bringing the stress term over from the RHS.

         # The gradient of the height of the free surface, h
         C_momentum = 0
         for dim in range(dimension):
            C_momentum += -g_magnitude*inner(self.w[dim], grad(self.h)[dim])*dx
         F -= C_momentum

         # Quadratic drag term in the momentum equation
         if(self.options["have_bottom_drag"]):
            print "Adding bottom drag..."
            
            # Get the drag coefficient C_D.
            # FIXME: This should be moved outside of the time loop.
            expr = ScalarExpressionFromOptions(path = "/equation/momentum_equation/drag_term/scalar_field::DragCoefficient/value", t=t)
            C_D = Function(self.W.sub(dimension)).interpolate(Expression(expr.code[0]))
            
            D_momentum = 0
            
            # Magnitude of the velocity field
            magnitude = 0
            for dim in range(dimension):
               magnitude += dot(self.u[dim], self.u[dim])
            magnitude = sqrt(magnitude)
            
            # Form the drag term
            for dim in range(dimension):
               D_momentum += -inner(self.w[dim], (C_D*magnitude/H)*self.u[dim])*dx
            F -= D_momentum

         # The mass term in the shallow water continuity equation 
         # (i.e. an advection equation for the free surface height, h)
         M_continuity = (1.0/dt)*(inner(self.v, self.h) - inner(self.v, self.h_old))*dx
         F += M_continuity

         # Divergence term in the shallow water continuity equation
         if(self.options["integrate_continuity_equation_by_parts"]):
            Ct_continuity = 0

            for dim in range(dimension):
               Ct_continuity += - H*inner(self.u[dim], grad(self.v)[dim])*dx
                                #+ inner(jump(v, n)[dim], avg(u[dim]))*dS
                              
            # Add in the surface integrals, but check to see if any boundary conditions need to be applied weakly here.
            boundary_markers = set(self.mesh.exterior_facets.markers)
            for marker in boundary_markers:
               marker = int(marker) # ds() will not accept markers of type 'numpy.int32', so convert it to type 'int' here.
               
               bc_type = None
               for i in range(0, libspud.option_count("/core_fields/vector_field::Velocity/boundary_conditions")):
                  bc_path = "/core_fields/vector_field::Velocity/boundary_conditions[%d]" % i
                  if(not (marker in libspud.get_option(bc_path + "/surface_ids"))):
                     # This BC is not associated with this marker, so skip it.
                     continue
                     
                  # Determine the BC type.
                  if(libspud.have_option(bc_path + "/type::no_normal_flow")):
                     bc_type = "no_normal_flow"
                  elif(libspud.have_option(bc_path + "/type::dirichlet")):
                     if(libspud.have_option(bc_path + "/type::dirichlet/apply_weakly")):
                        bc_type = "weak_dirichlet"
                     else:
                        bc_type = "dirichlet"
                  elif(libspud.have_option(bc_path + "/type::flather")):
                     bc_type = "flather"
                     
                  # Apply the boundary condition...
                  if(bc_type == "flather"):
                     print "Applying flather BC to surface ID %d..." % marker
                     
                     # The known exterior value for the Velocity.
                     expr = VectorExpressionFromOptions(path = (bc_path + "/type::flather/exterior_velocity"), t=t)
                     u_ext = [Function(self.W.sub(dim)).interpolate(Expression(expr.code[dim])) for dim in range(dimension)]
                     for dim in range(dimension):
                        Ct_continuity += H*inner(u_ext[dim], self.n[dim]) * self.v * ds(int(marker))
                     
                     # The known exterior value for the FreeSurfacePerturbation.
                     expr = ScalarExpressionFromOptions(path = (bc_path + "/type::flather/exterior_free_surface_perturbation"), t=t)
                     h_ext = Function(self.W.sub(dimension)).interpolate(Expression(expr.code[0]))
                     Ct_continuity += H*sqrt(g_magnitude/H)*inner(self.h - h_ext, self.v) * ds(int(marker))
                     
                  elif(bc_type == "weak_dirichlet"):
                     print "Applying weak Dirichlet BC to surface ID %d..." % marker
                     expr = VectorExpressionFromOptions(path = (bc_path + "/type::dirichlet"), t=t)
                     u_bdy = [Function(self.W.sub(dim)).interpolate(Expression(expr.code[dim])) for dim in range(dimension)]
                     for dim in range(dimension):
                        Ct_continuity += H*inner(u_bdy[dim], self.n[dim]) * self.v * ds(int(marker))
                  elif(bc_type == "dirichlet"):
                     # Add in the surface integral as it is here. The BC will be applied strongly later using a DirichletBC object.
                     for dim in range(dimension):
                        Ct_continuity += H*inner(self.u[dim], self.n[dim]) * self.v * ds(int(marker))
                  elif(bc_type == "no_normal_flow"):
                     print "Applying no normal flow BC to surface ID %d..." % marker
                     # Do nothing here since dot(u, n) is zero.
                  else:
                     print "Unknown boundary condition type!"
                     sys.exit(0)
                     
               # If no boundary condition has been applied, include the surface integral as it is.
               if(bc_type is None):
                  for dim in range(0, dimension):
                     Ct_continuity += H*inner(self.u[dim], self.n[dim]) * self.v * ds(int(marker))

         else:
            divergence = 0
            for dim in range(dimension):
               divergence += grad(H*self.u[dim])[dim]
            Ct_continuity = inner(self.v, divergence)*dx
         F += Ct_continuity
            
         # Add in any source terms
         if(self.options["have_momentum_source"]):
            expr = VectorExpressionFromOptions(path = "/equations/momentum_equation/source_term/vector_field::Source/value", t=t)
            momentum_source = [Function(self.W.sub(dim)).interpolate(Expression(expr.code[dim])) for dim in range(dimension)]
            print "Adding momentum source..."
            for dim in range(dimension):
               F -= inner(self.w[dim], momentum_source[dim])*dx

         if(self.options["have_continuity_source"]):
            expr = ScalarExpressionFromOptions(path = "/equations/continuity_equation/source_term/vector_field::Source/value", t=t)
            continuity_source = Function(self.W.sub(dimension)).interpolate(Expression(expr.code[0]))
            print "Adding continuity source..."
            F -= inner(self.v, continuity_source)*dx
            
         # Add in any SU stabilisation
         if(self.options["have_su_stabilisation"]):
            print "Adding momentum SU stabilisation..."
            stabilisation = Stabilisation(self.mesh, P1, cellsize)
            u_temp = []
            for dim in range(dimension):
               u_temp.append(self.solution_old.split()[dim])
            F += stabilisation.streamline_upwind(self.w, self.u, u_temp, self.options["viscosity"])

         # Get all the Dirichlet boundary conditions for the Velocity field
         bcs = []
         for i in range(0, libspud.option_count("/core_fields/vector_field::Velocity/boundary_condition")):
            if(libspud.have_option("/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet" % i) and
               not libspud.have_option("/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet/apply_weakly" % i)):
               # If it's not a weak BC, then it must be a strong one.
               print "Applying Velocity BC #%d" % i
               expr = VectorExpressionFromOptions(path = ("/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet" % i), t=t)
               # Surface IDs on the domain boundary
               surface_ids = libspud.get_option("/core_fields/vector_field::Velocity/boundary_condition[%d]/surface_ids" % i)
               for dim in range(dimension):
                  bc = DirichletBC(self.W.sub(dim), Expression(expr.code[dim]), surface_ids)
                  bcs.append(bc)

         # Get all the Dirichlet boundary conditions for the FreeSurfacePerturbation field
         for i in range(0, libspud.option_count("/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition/type::dirichlet")):
            if(libspud.have_option("/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/type::dirichlet" % i) and
               not(libspud.have_option("/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/type::dirichlet/apply_weakly" % i))):
               # If it's not a weak BC, then it must be a strong one.
               print "Applying FreeSurfacePerturbation BC #%d" % i
               expr = ScalarExpressionFromOptions(path = ("/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/type::dirichlet" % i), t=t)
               # Surface IDs on the domain boundary
               surface_ids = libspud.get_option("/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/surface_ids" % i)
               bc = DirichletBC(self.W.sub(dimension), Expression(expr.code[0]), surface_ids)
               bcs.append(bc)
            
         # Split up the form into the LHS and RHS (i.e. the bilinear and linear forms)
         #a = lhs(F)
         #L = rhs(F)
                   
         # Solve the system of equations!
         #solve(A, self.solution, b, bcs=bcs, solver_parameters={'ksp_monitor':True, 'ksp_view':True, 'pc_view':True})
         solve(F == 0, self.solution, bcs=bcs, solver_parameters={'ksp_monitor': True, 
                                                                  'ksp_view': False, 
                                                                  'pc_view': False, 
                                                                  'ksp_type': 'gmres', 
                                                                  'pc_type': 'jacobi',
                                                                  'ksp_rtol': 1.0e-7,
                                                                  'snes_type': 'ksponly'})
            
         # Write the solution to file.
         print "Writing data to file..."
         for dim in range(dimension):
            projected = project(self.solution.split()[dim], self.W.sub(dimension))
            self.output_function[dim].assign(projected)
            self.output_file[dim] << self.output_function[dim]
         self.output_function[dimension].assign(self.solution.split()[dimension])
         self.output_file[dimension] << self.output_function[dimension]
         
         # Check whether a steady-state has been reached.
         # Take the maximum difference across all processes.
         global_max_difference_h = max(abs(self.solution.split()[dimension].vector().gather() - self.solution_old.split()[dimension].vector().gather()))
         global_max_difference_u = [max(abs(self.solution.split()[dim].vector().gather() - self.solution_old.split()[dim].vector().gather())) for dim in range(0, dimension)]
         # If the difference is less than a set tolerance, then break out of the time-stepping loop.
         if(global_max_difference_h <= self.options["steady_state_tolerance"] and (numpy.array(global_max_difference_u) <= self.options["steady_state_tolerance"]).all()):
            print "Steady-state attained. Exiting the time-stepping loop..."
            break

         # Write detector values to file
         if(self.options["have_detectors"]):
            self.detectors.write(self.options["simulation_name"], t, dt)
            
         # Move to next time step    
         self.solution_old.assign(self.solution)    
         t += dt
         print "Moving to next time level..."      
      
      print "Out of the time-stepping loop."
   
      # Any final steps (e.g. closing files)
      if(self.options["have_detectors"]):
         self.detectors.finalise()
  
      return

if(__name__ == "__main__"):

   try:
      path = sys.argv[1]
   except IndexError:
      print "Please provide the path to the simulation configuration file."
      sys.exit(1)
      
   # Set up a shallow water simulation.
   sw = ShallowWater(path)
      
   # Solve the shallow water equations!
   sw.run()
   

