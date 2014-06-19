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

class ExpressionFromOptions(Expression):
   def __init__(self, path, t):
      self.path = path
      if(libspud.have_option(path + "/constant")):
         self.source_type = "constant"
         self.source_value = libspud.get_option(path + "/constant")
         Expression.__init__(self, code=self.source_value, t=t)  
      elif(libspud.have_option(path + "/cpp")):
         self.source_type = "cpp"
         self.source_value = libspud.get_option(path + "/cpp")
         exec self.source_value
         Expression.__init__(self, code=val(t), t=t)
      else:
         print "No value specified"
         sys.exit(1)

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
            if(name == "VelocityFunctionSpace"): # FIXME: Include a "scalar/vector/tensor" type for the function space in the schema.
               self.function_spaces[name] = VectorFunctionSpace(self.mesh, "CG", degree)
            else:
               self.function_spaces[name] = FunctionSpace(self.mesh, "CG", degree)
         elif(family == "Discontinuous Lagrange"):
            if(name == "VelocityFunctionSpace"):
               self.function_spaces[name] = VectorFunctionSpace(self.mesh, "DG", degree)
            else:
               self.function_spaces[name] = FunctionSpace(self.mesh, "DG", degree)
         else:
            print "Unknown element family: %s." % family
            sys.exit(1)
            
      # Define the coordinate function space as P1.
      self.function_spaces["CoordinateFunctionSpace"] = FunctionSpace(self.mesh, "CG", 1)
                
      # Define the mixed function space
      U = self.function_spaces["VelocityFunctionSpace"]
      H = self.function_spaces["FreeSurfaceFunctionSpace"]
      self.W = MixedFunctionSpace([U, H])
     
      # The solution field defined on the mixed function space
      self.solution = Function(self.W)
      # These are like the TrialFunctions, but are just regular Functions here because we want to solve a non-linear problem
      functions = split(self.solution)
      self.u = functions[0]; self.h = functions[1]

      # Get the test functions
      test_functions = TestFunctions(self.W)
      self.w = test_functions[0]; self.v = test_functions[1]

      # Normal vector to each element facet
      self.n = FacetNormal(self.mesh)
      
      # Set up initial conditions
      # FIXME: Subclassing of the Expression class needs to be DOLFIN compatible. The current method used here is a hack.
      h_initial = ExpressionFromOptions(path = "/system/core_fields/scalar_field::FreeSurfacePerturbation/initial_condition", t=self.options["t"])
      expr = ExpressionFromOptions(path = "/system/core_fields/vector_field::Velocity/initial_condition", t=self.options["t"])
      u_initial = [expr.code[dim] for dim in range(dimension)]

      # Define the compulsory shallow water fields
      codes = tuple(u_initial) + (h_initial.code[0],)
      self.solution_old = Function(self.W).interpolate(Expression(codes))
      functions_old = split(self.solution_old)
      self.u_old = functions_old[0]; self.h_old = functions_old[1]
      
      # The solution should first hold the initial condition.
      self.solution.assign(self.solution_old)
      
      # Mean free surface height
      self.h_mean = Function(self.W.sub(1))
      self.h_mean.interpolate(ExpressionFromOptions(path = "/system/core_fields/scalar_field::FreeSurfaceMean/value", t=self.options["t"]))

      # Set up the functions used to write fields to file.
      self.output_function = [Function(self.W.sub(0), name="Velocity"), Function(self.W.sub(1), name="FreeSurfacePerturbation")]
      
      # Set up the output stream
      self.output_file = [File("%s_Velocity.pvd" % self.options["simulation_name"]), File("%s_FreeSurfacePerturbation.pvd" % self.options["simulation_name"])]
      
      # Write initial conditions to file
      self.output_function[0].assign(self.solution_old.split()[0])
      self.output_file[0] << self.output_function[0]
      self.output_function[1].assign(self.solution_old.split()[1])
      self.output_file[1] << self.output_function[1]
    
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
      if(libspud.have_option("/system/equations/momentum_equation/mass_term/exclude_mass_term")):
         self.options["have_momentum_mass"] = False
      else:
         self.options["have_momentum_mass"] = True
         
      if(libspud.have_option("/system/equations/momentum_equation/advection_term/exclude_advection_term")):
         self.options["have_momentum_advection"] = False
      else:
         self.options["have_momentum_advection"] = True
         
      self.options["have_momentum_stress"] = libspud.have_option("/system/equations/momentum_equation/stress_term")
         
      # Source terms for the momentum and continuity equations
      self.options["have_momentum_source"] = libspud.have_option("/system/equations/momentum_equation/source_term")
      self.options["have_continuity_source"] = libspud.have_option("/system/equations/continuity_equation/source_term")

      # Check for any SU stabilisation
      self.options["have_su_stabilisation"] = libspud.have_option("/system/equations/momentum_equation/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind_stabilisation")
      
      # Turbulence parameterisation
      self.options["have_turbulence_parameterisation"] = libspud.have_option("/system/equations/momentum_equation/turbulence_parameterisation")

      self.options["have_bottom_drag"] = libspud.have_option("/system/equations/momentum_equation/drag_term")
         
      self.options["integrate_continuity_equation_by_parts"] = libspud.have_option("/system/equations/continuity_equation/integrate_by_parts")
      self.options["integrate_advection_term_by_parts"] = libspud.have_option("/system/equations/momentum_equation/advection_term/integrate_by_parts")

      self.options["have_detectors"] = libspud.have_option("/io/detectors/")
      
      return
      
   def run(self):
      """ Execute the time-stepping loop. """
   
      # Time-stepping parameters and constants
      T = self.options["T"]
      t = self.options["t"]
      dt = self.options["dt"]
      dimension = self.options["dimension"]
      g_magnitude = self.options["g_magnitude"]
      
      # The total height of the free surface.
      H = self.h_mean + self.h
      
      # Coordinate mesh
      P1 = self.function_spaces["CoordinateFunctionSpace"]
      cellsize = CellSize(self.mesh)

      # The collection of all the individual terms in their weak form.
      F = 0

      # Mass term
      if(self.options["have_momentum_mass"]):
         print "Adding mass term..."
         M_momentum = (1.0/dt)*(inner(self.w, self.u) - inner(self.w, self.u_old))*dx
         F += M_momentum
      
      # Advection term
      if(self.options["have_momentum_advection"]):
         print "Adding advection term..."

         if(self.options["integrate_advection_term_by_parts"]):
            outflow = (dot(self.u, self.n) + abs(dot(self.u, self.n)))/2.0

            A_momentum = -inner(dot(self.u, grad(self.w)), self.u)*dx - inner(dot(self.u, grad(self.u)), self.w)*dx
            A_momentum += inner(self.w, outflow*self.u)*ds
            A_momentum += dot(outflow('+')*self.u('+') - outflow('-')*self.u('-'), jump(self.w))*dS

         else:
            A_momentum = inner(dot(grad(self.u), self.u), self.w)*dx
         F += A_momentum
         
      # Viscous stress term. Note that the viscosity is kinematic (not dynamic).
      if(self.options["have_momentum_stress"]):
         print "Adding stress term..."
         
         # Background viscosity
         viscosity = Constant(libspud.get_option("/system/equations/momentum_equation/stress_term/scalar_field::Viscosity/value/constant"))
         
         # Eddy viscosity
         if(self.options["have_turbulence_parameterisation"]):
            base_option_path = "/system/equations/momentum_equation/turbulence_parameterisation"
            # Large eddy simulation (LES)
            if(libspud.have_option(base_option_path + "/les")):
               les = LES(self.mesh, self.W.sub(1))
               density = Function(self.W.sub(1)).interpolate(Constant("1.0")) # We divide through by density in the momentum equation, so just set this to 1.0 for now.
               smagorinsky_coefficient = Function(self.W.sub(1)).interpolate(Constant(libspud.get_option(base_option_path + "/les/smagorinsky/smagorinsky_coefficient")))
               filter_width = Function(self.W.sub(1)).interpolate(Constant(libspud.get_option(base_option_path + "/les/smagorinsky/filter_width"))) # FIXME: Remove this when CellSize is supported in Firedrake.
               eddy_viscosity = les.eddy_viscosity(self.u, density, smagorinsky_coefficient, filter_width)
            # Add on eddy viscosity
            viscosity += eddy_viscosity
            
         # Stress tensor: tau = grad(u) + transpose(grad(u)) - (2/3)*div(u)
         dg = False
         if(not dg):
            # Perform a double dot product of the stress tensor and grad(w).
            K_momentum = -viscosity*inner(grad(self.u) + grad(self.u).T, grad(self.w))*dx
            #K_momentum += viscosity*(2.0/3.0)*inner(div(self.u)*Identity(dimension), grad(self.w))*dx
         else:
            # Interior penalty method
            cellsize = Constant(0.2) #CellSize(self.mesh)
            alpha = Constant(5.0) # Penalty parameter.
            
            K_momentum = -viscosity('+')*inner(grad(self.u), grad(self.w))*dx
            for dim in range(len(self.u)):
               K_momentum += -viscosity('+')*(alpha('+')/cellsize('+'))*dot(jump(self.w[dim], self.n), jump(self.u[dim], self.n))*dS
               K_momentum += viscosity('+')*dot(avg(grad(self.w[dim])), jump(self.u[dim], self.n))*dS + viscosity('+')*dot(jump(self.w[dim], self.n), avg(grad(self.u[dim])))*dS

         F -= K_momentum # Negative sign here because we are bringing the stress term over from the RHS.

      # The gradient of the height of the free surface, h
      C_momentum = -g_magnitude*inner(self.w, grad(self.h))*dx
      F -= C_momentum

      # Quadratic drag term in the momentum equation
      if(self.options["have_bottom_drag"]):
         print "Adding bottom drag..."
         
         # Get the drag coefficient C_D.
         C_D = Function(self.W.sub(1)).interpolate(ExpressionFromOptions(path="/system/equations/momentum_equation/drag_term/scalar_field::DragCoefficient/value", t=t))
         
         # Magnitude of the velocity field
         magnitude = sqrt(dot(self.u, self.u))
         
         # Form the drag term
         D_momentum = -inner(self.w, (C_D*magnitude/H)*self.u)*dx
         F -= D_momentum

      # The mass term in the shallow water continuity equation 
      # (i.e. an advection equation for the free surface height, h)
      M_continuity = (1.0/dt)*(inner(self.v, self.h) - inner(self.v, self.h_old))*dx
      F += M_continuity

      # Append any Expression objects for weak BCs here.
      weak_bc_expressions = []
      
      # Divergence term in the shallow water continuity equation
      if(self.options["integrate_continuity_equation_by_parts"]):

         Ct_continuity = - H*inner(self.u, grad(self.v))*dx \
                         + inner(jump(self.v, self.n), avg(H*self.u))*dS
                           
         # Add in the surface integrals, but check to see if any boundary conditions need to be applied weakly here.
         boundary_markers = self.mesh.exterior_facets.unique_markers
         for marker in boundary_markers:
            marker = int(marker) # ds() will not accept markers of type 'numpy.int32', so convert it to type 'int' here.
            
            bc_type = None
            for i in range(0, libspud.option_count("/system/core_fields/vector_field::Velocity/boundary_condition")):
               bc_path = "/system/core_fields/vector_field::Velocity/boundary_condition[%d]" % i
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
                  u_ext = ExpressionFromOptions(path = (bc_path + "/type::flather/exterior_velocity"), t=t)
                  #expr = Constant([x_val, y_val]) 
                  Ct_continuity += H*inner(Function(self.W.sub(0)).interpolate(u_ext), self.n)*self.v*ds(int(marker))
                  
                  # The known exterior value for the FreeSurfacePerturbation.
                  h_ext = ExpressionFromOptions(path = (bc_path + "/type::flather/exterior_free_surface_perturbation"), t=t)
                  Ct_continuity += H*sqrt(g_magnitude/H)*(self.h - Function(self.W.sub(1)).interpolate(h_ext))*self.v*ds(int(marker))
                  
                  weak_bc_expressions.append(u_ext)
                  weak_bc_expressions.append(h_ext)
                  
               elif(bc_type == "weak_dirichlet"):
                  print "Applying weak Dirichlet BC to surface ID %d..." % marker
                  u_bdy = ExpressionFromOptions(path = (bc_path + "/type::dirichlet"), t=t)
                  Ct_continuity += H*(dot(Function(self.W.sub(0)).interpolate(u_bdy), self.n))*self.v*ds(int(marker))
                  
                  weak_bc_expressions.append(u_bdy)
                  
               elif(bc_type == "dirichlet"):
                  # Add in the surface integral as it is here. The BC will be applied strongly later using a DirichletBC object.
                  Ct_continuity += H * inner(self.u, self.n) * self.v * ds(int(marker))
               elif(bc_type == "no_normal_flow"):
                  print "Applying no normal flow BC to surface ID %d..." % marker
                  # Do nothing here since dot(u, n) is zero.
               else:
                  print "Unknown boundary condition type!"
                  sys.exit(0)
                  
            # If no boundary condition has been applied, include the surface integral as it is.
            if(bc_type is None):
               Ct_continuity += H * inner(self.u, self.n) * self.v * ds(int(marker))

      else:
         Ct_continuity = inner(self.v, div(H*self.u))*dx
      F += Ct_continuity

      # Add in any source terms
      if(self.options["have_momentum_source"]):
         print "Adding momentum source..."
         momentum_source = ExpressionFromOptions(path = "/system/equations/momentum_equation/source_term/vector_field::Source/value", t=t)
         F -= inner(self.w, Function(self.W.sub(0)).interpolate(momentum_source))*dx

      if(self.options["have_continuity_source"]):
         print "Adding continuity source..."
         continuity_source = ExpressionFromOptions(path = "/system/equations/continuity_equation/source_term/scalar_field::Source/value", t=t)
         F -= inner(self.v, Function(self.W.sub(1)).interpolate(continuity_source))*dx
         

      # Add in any SU stabilisation
      if(self.options["have_su_stabilisation"]):
         print "Adding momentum SU stabilisation..."
         stabilisation = Stabilisation(self.mesh, P1, cellsize)
         
         magnitude = fields_calculations.magnitude_vector(self.mesh, self.solution_old.split()[0], P1)

         # Bound the values for the magnitude below by 1.0e-9 for numerical stability reasons.
         u_nodes = magnitude.vector()
         near_zero = numpy.array([1.0e-9 for i in range(len(u_nodes))])
         u_nodes.set_local(numpy.maximum(u_nodes.array(), near_zero))

         diffusivity = ExpressionFromOptions(path = "/system/equations/momentum_equation/stress_term/scalar_field::Viscosity/value", t=self.options["t"])
         diffusivity = Function(self.W.sub(1)).interpolate(diffusivity) # Background viscosity
         grid_pe = fields_calculations.grid_peclet_number(self.mesh, diffusivity, magnitude, P1, cellsize)
   
         # Bound the values for grid_pe below by 1.0e-9 for numerical stability reasons. 
         grid_pe_nodes = grid_pe.vector()
         values = numpy.array([1.0e-9 for i in range(len(grid_pe_nodes))])
         grid_pe_nodes.set_local(numpy.maximum(grid_pe_nodes.array(), values))

         F += stabilisation.streamline_upwind(self.w, self.u, magnitude, grid_pe)

      # Get all the Dirichlet boundary conditions for the Velocity field
      bcs = []
      bc_expressions = []
      for i in range(0, libspud.option_count("/system/core_fields/vector_field::Velocity/boundary_condition")):
         if(libspud.have_option("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet" % i) and
            not libspud.have_option("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet/apply_weakly" % i)):
            # If it's not a weak BC, then it must be a strong one.
            print "Applying Velocity BC #%d" % i
            expr = ExpressionFromOptions(path = ("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet" % i), t=t)
            # Surface IDs on the domain boundary
            surface_ids = libspud.get_option("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/surface_ids" % i)
            bc = DirichletBC(self.W.sub(0), expr, surface_ids, method="geometric")
            bcs.append(bc)
            bc_expressions.append(expr)

      # Get all the Dirichlet boundary conditions for the FreeSurfacePerturbation field
      for i in range(0, libspud.option_count("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition/type::dirichlet")):
         if(libspud.have_option("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/type::dirichlet" % i) and
            not(libspud.have_option("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/type::dirichlet/apply_weakly" % i))):
            # If it's not a weak BC, then it must be a strong one.
            print "Applying FreeSurfacePerturbation BC #%d" % i
            expr = ExpressionFromOptions(path = ("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/type::dirichlet" % i), t=t)
            # Surface IDs on the domain boundary
            surface_ids = libspud.get_option("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/surface_ids" % i)
            bc = DirichletBC(self.W.sub(1), expr, surface_ids, method="geometric")
            bcs.append(bc)
            bc_expressions.append(expr)

      # Construct the solver objects
      problem = NonlinearVariationalProblem(F, self.solution, bcs=bcs)
      solver = NonlinearVariationalSolver(problem, solver_parameters={'ksp_monitor': True, 
                                                                  'ksp_view': False, 
                                                                  'pc_view': False, 
                                                                  'pc_type': 'fieldsplit',
                                                                  'pc_fieldsplit_type': 'schur',
                                                                  'ksp_type': 'gmres',
                                                                  'pc_fieldsplit_schur_fact_type': 'FULL',
                                                                  'fieldsplit_0_ksp_type': 'preonly',
                                                                  'fieldsplit_1_ksp_type': 'preonly',
                                                                  'ksp_rtol': 1.0e-7,
                                                                  'snes_type': 'ksponly'})
      
      t += dt

      # The time-stepping loop
      EPSILON = 1.0e-14
      while t <= T + EPSILON: # A small value EPSILON is added here in case of round-off error.
         print "\nt = %g" % t

         # Update any time-dependent Functions and Expressions.
         
         # Re-compute the velocity magnitude and grid Peclet number fields.
         if(self.options["have_su_stabilisation"]):
            magnitude.assign(fields_calculations.magnitude_vector(self.mesh, self.solution_old.split()[0], P1))

            # Bound the values for the magnitude below by 1.0e-9 for numerical stability reasons.
            u_nodes = magnitude.vector()
            near_zero = numpy.array([1.0e-9 for i in range(len(u_nodes))])
            u_nodes.set_local(numpy.maximum(u_nodes.array(), near_zero))
            
            grid_pe.assign(fields_calculations.grid_peclet_number(self.mesh, diffusivity, magnitude, P1, cellsize))
      
            # Bound the values for grid_pe below by 1.0e-9 for numerical stability reasons. 
            grid_pe_nodes = grid_pe.vector()
            values = numpy.array([1.0e-9 for i in range(len(grid_pe_nodes))])
            grid_pe_nodes.set_local(numpy.maximum(grid_pe_nodes.array(), values))

         # Time-dependent source terms
         if(self.options["have_momentum_source"]):
            momentum_source.t = t
         if(self.options["have_continuity_source"]):
            continuity_source.t = t
            
         # Update any time-varying DirichletBC objects.
         for expr in bc_expressions:
            expr.t = t
         for expr in weak_bc_expressions:
            expr.t = t
                   
         # Solve the system of equations!
         solver.solve()
            
         # Write the solution to file.
         print "Writing data to file..."
         self.output_function[0].assign(self.solution.split()[0])
         self.output_file[0] << self.output_function[0]
         self.output_function[1].assign(self.solution.split()[1])
         self.output_file[1] << self.output_function[1]
         
         # Check whether a steady-state has been reached.
         # Take the maximum difference across all processes.
         global_max_difference_h = max(abs(self.solution.split()[1].vector().gather() - self.solution_old.split()[1].vector().gather()))
         global_max_difference_u = max(abs(self.solution.split()[0].vector().gather() - self.solution_old.split()[0].vector().gather()))
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
   

