#    Copyright (C) 2013 Imperial College London.

#    This file is part of Firedrake-Fluids.
#
#    Firedrake-Fluids is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Firedrake-Fluids is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Firedrake-Fluids.  If not, see <http://www.gnu.org/licenses/>.

import sys, os
from firedrake_fluids import LOG

# PETSc environment variables
try:
   if(os.environ["PETSC_OPTIONS"] == ""):
      os.environ["PETSC_OPTIONS"] = "-log_summary"
   else:
      os.environ["PETSC_OPTIONS"] = os.environ["PETSC_OPTIONS"] + " -log_summary"
except KeyError:
   # Environment variable does not exist, so let's set it now.
   os.environ["PETSC_OPTIONS"] = "-log_summary"
LOG.debug("Environment variable PETSC_OPTIONS set to: %s" % (os.environ["PETSC_OPTIONS"]))

import numpy
import mpi4py

import libspud
LOG.debug("libspud successfully imported")
from pyop2 import *
LOG.debug("PyOP2 successfully imported")
from firedrake import *
LOG.debug("Firedrake successfully imported")

# COFFEE, PyOP2 and FFC parameters
op2.init(lazy_evaluation=False)
parameters['form_compiler']['quadrature_degree'] = 4
parameters["coffee"]["O2"] = False # FIXME: Remove this one this issue has been fixed: https://github.com/firedrakeproject/firedrake/issues/425

# Firedrake-Fluids modules
from firedrake_fluids.utils import *
from firedrake_fluids.fields_calculations import *
from firedrake_fluids.stabilisation import Stabilisation
from firedrake_fluids.les import LES
from firedrake_fluids.expression import ExpressionFromOptions
from firedrake_fluids.metadata import *
from firedrake_fluids.diagnostics import Diagnostics
LOG.debug("Firedrake-Fluids sub-modules successfully imported.")

class ShallowWater:
   r""" A class for setting up and running a non-linear shallow water simulation.
   
   At its core, this involves solving a coupled system of equations for two primitive variables :math:`\mathbf{u}` (the velocity) and :math:`h` (the perturbation of the free surface).
   
   The momentum equation is given by
   
   .. math:: \frac{\partial \mathbf{u}}{\partial t} + \mathbf{u}\cdot\nabla\mathbf{u} = -g\nabla h + \nabla\cdot\mathbb{T} - C_D\frac{||\mathbf{u}||\mathbf{u}}{(H + h)},

   where :math:`g` is the acceleration due to gravity, :math:`\mathbf{u}`
   is the velocity, and :math:`C_D` is the non-dimensional drag
   coefficient. The stress tensor :math:`\mathbb{T}` is given by

   .. math:: \mathbb{T} = \nu\left(\nabla\mathbf{u} + \nabla\mathbf{u}^{\mathrm{T}}\right) - \frac{2}{3}\nu\left(\nabla\cdot\mathbf{u}\right)\mathbb{I},

   where :math:`\nu` is the isotropic kinematic viscosity, and
   :math:`\mathbb{I}` is the identity tensor.

   The continuity equation is given by

   .. math:: \frac{\partial h}{\partial t} + \nabla\cdot\left(\left(H + h\right)\mathbf{u}\right) = 0.   
   """
   
   def __init__(self, path, checkpoint=None):
      """ Initialise a new shallow water simulation, using an options file.
      
      :param str path: The path to the simulation's configuration/options file.
      :param str checkpoint: The path to a checkpoint file.
      :returns: None
      """
   
      LOG.info("Initialising simulation...")
      
      # Remove any stored options.
      libspud.clear_options()

      # Load the options from the options tree.
      libspud.load_options(path)
     
      # Populate the options dictionary
      self.populate_options()
        
      # Read in the input mesh (or construct one)
      self.mesh = self.get_mesh(path)

      # Create a dictionary containing all the function spaces
      LOG.info("Creating function spaces...")
      self.function_spaces = {}
                      
      # The Velocity field's function space
      name = "VelocityFunctionSpace"
      path = "/function_spaces/function_space::%s" % name
      family = libspud.get_option(path+"/family")
      degree = libspud.get_option(path+"/degree")
      try:
         if(family == "Continuous Lagrange"):
            self.function_spaces[name] = VectorFunctionSpace(self.mesh, "CG", degree)
         elif(family == "Discontinuous Lagrange"):
            self.function_spaces[name] = VectorFunctionSpace(self.mesh, "DG", degree)
         else:
            raise ValueError("Unknown element family: %s." % family)
      except ValueError as e:      
         LOG.exception(e)
         sys.exit()
      LOG.debug("Created a new %s function space of degree %d for Velocity" % (family, degree))
      
      # The FreeSurfacePerturbation field's function space
      name = "FreeSurfaceFunctionSpace"
      path = "/function_spaces/function_space::%s" % name
      family = libspud.get_option(path+"/family")
      degree = libspud.get_option(path+"/degree")
      try:
         if(family == "Continuous Lagrange"):
            self.function_spaces[name] = FunctionSpace(self.mesh, "CG", degree)
         elif(family == "Discontinuous Lagrange"):
            self.function_spaces[name] = FunctionSpace(self.mesh, "DG", degree)
         else:
            raise ValueError("Unknown element family: %s." % family)
      except ValueError as e:
         LOG.exception(e)
         sys.exit()
      LOG.debug("Created a new %s function space of degree %d for FreeSurfacePerturbation" % (family, degree))
      
      # Get the function spaces
      U = self.function_spaces["VelocityFunctionSpace"]
      H = self.function_spaces["FreeSurfaceFunctionSpace"]

      # Test functions
      self.w = TestFunction(U)
      self.v = TestFunction(H)
      LOG.info("Test functions created.")
      
      # Trial functions (i.e. the solutions for time n+1)
      self.u = TrialFunction(U)
      self.h = TrialFunction(H)
      LOG.info("Trial functions created.")

      # Normal vector to each element facet
      self.n = FacetNormal(self.mesh)

      # Set up initial conditions
      LOG.info("Setting initial conditions...")
      h_initial = ExpressionFromOptions(path = "/system/core_fields/scalar_field::FreeSurfacePerturbation/initial_condition").get_expression()
      u_initial = ExpressionFromOptions(path = "/system/core_fields/vector_field::Velocity/initial_condition").get_expression()
      
      # The solution at time n.
      self.u0 = Function(U).interpolate(u_initial)
      self.h0 = Function(H).interpolate(h_initial)
      
      # The solution at time n-1.
      self.u00 = Function(U)
      self.u00.assign(self.u0)
      
      # Load initial conditions from the specified checkpoint file if desired.
      if(checkpoint is not None):
         self.solution_old.dat.load(checkpoint)
         LOG.debug("Loaded initial condition from checkpoint file.")
      
      # Mean free surface height
      self.h_mean = Function(H)
      self.h_mean.interpolate(ExpressionFromOptions(path = "/system/core_fields/scalar_field::FreeSurfaceMean/value").get_expression())

      # Set up the functions used to write fields to file.
      self.output_functions = {}
      self.output_functions["Velocity"] = Function(U, name="Velocity")
      self.output_functions["FreeSurfacePerturbation"] = Function(H, name="FreeSurfacePerturbation")
      
      # Set up the output stream
      LOG.info("Initialising output file streams...")
      self.output_files = {}
      for field in self.output_functions.keys():
         self.output_files[field] = File("%s_%s.pvd" % (self.options["simulation_name"], field))
      
      # Write initial conditions to file
      LOG.info("Writing initial conditions to file...")
      self.output_functions["Velocity"].assign(self.u0)
      self.output_files["Velocity"] << self.output_functions["Velocity"]
      self.output_functions["FreeSurfacePerturbation"].assign(self.h0)
      self.output_files["FreeSurfacePerturbation"] << self.output_functions["FreeSurfacePerturbation"]
      
      return
      
   def populate_options(self):
      """ Add simulation options related to the shallow water model to a dictionary object.
      
      :returns: None
      """
      # A dictionary storing all the options
      self.options = {}
      
      LOG.info("Populating options dictionary...")
      
      self.options["simulation_name"] = libspud.get_option("/simulation_name")
      self.options["dimension"] = libspud.get_option("/geometry/dimension")
   
      # Time-stepping parameters
      self.options["T"] = libspud.get_option("/timestepping/finish_time")
      self.options["t"] = libspud.get_option("/timestepping/current_time")
      self.options["dt"] = libspud.get_option("/timestepping/timestep")
      self.options["theta"] = libspud.get_option("/timestepping/theta")
      if(libspud.have_option("/timestepping/steady_state")):
         self.options["steady_state_tolerance"] = libspud.get_option("/timestepping/steady_state/tolerance")
      else:
         self.options["steady_state_tolerance"] = -1000
         
      # I/O parameters
      if(libspud.have_option("/io/dump_period")):
         self.options["dump_period"] = libspud.get_option("/io/dump_period")
      else:
         self.options["dump_period"] = None
     
      if(libspud.have_option("/io/checkpoint")):
         self.options["checkpoint_period"] = libspud.get_option("/io/checkpoint/dump_period")
      else:
         self.options["checkpoint_period"] = None
           
      # Physical parameters
      self.options["g_magnitude"] = libspud.get_option("/physical_parameters/gravity/magnitude")
      
      # Enable/disable terms in the shallow water equations
      self.options["have_momentum_mass"] = (not libspud.have_option("/system/equations/momentum_equation/mass_term/exclude_mass_term"))
      self.options["have_momentum_advection"] = (not libspud.have_option("/system/equations/momentum_equation/advection_term/exclude_advection_term"))
      self.options["have_momentum_stress"] = libspud.have_option("/system/equations/momentum_equation/stress_term")
         
      # Source terms for the momentum and continuity equations
      self.options["have_momentum_source"] = libspud.have_option("/system/equations/momentum_equation/source_term")
      self.options["have_continuity_source"] = libspud.have_option("/system/equations/continuity_equation/source_term")

      # Check for any SU stabilisation
      self.options["have_su_stabilisation"] = libspud.have_option("/system/equations/momentum_equation/spatial_discretisation/continuous_galerkin/stabilisation/streamline_upwind")
      
      # Turbulence parameterisation
      self.options["have_turbulence_parameterisation"] = libspud.have_option("/system/equations/momentum_equation/turbulence_parameterisation")

      # Drag parameterisation
      self.options["have_drag"] = libspud.have_option("/system/equations/momentum_equation/drag_term")
         
      # Integration by parts
      self.options["integrate_continuity_equation_by_parts"] = libspud.have_option("/system/equations/continuity_equation/integrate_by_parts")
      self.options["integrate_advection_term_by_parts"] = libspud.have_option("/system/equations/momentum_equation/advection_term/integrate_by_parts")
      
      return

      
   def compute_diagnostics(self):
      diagnostic_field_count = libspud.option_count("/system/diagnostic_fields/diagnostic")
      if(diagnostic_field_count == 0):
         return
         
      LOG.info("Computing diagnostic fields...")
      d = Diagnostics(self.mesh)
      for i in range(0, diagnostic_field_count):
         name = libspud.get_option("/system/diagnostic_fields/diagnostic[%d]/name" % i)
         try:
            if(name == "grid_reynolds_number"):
               viscosity = Function(self.W.sub(1)).interpolate(Expression(libspud.get_option("/system/equations/momentum_equation/stress_term/scalar_field::Viscosity/value/constant")))
               field = d.grid_reynolds_number(self.u, viscosity)
            elif(name == "courant_number"):
               field = d.courant_number(self.u, self.options["dt"])
            else:
               raise ValueError("Unknown diagnostic field: %s" % name)
         except ValueError as e:
            LOG.exception(name)
            sys.exit()
            
         LOG.info("Diagnostic results for: %s" % name)
         LOG.info("Maximum value: %f" % max(field.vector()))
         LOG.info("Maximum value: %f" % min(field.vector()))
         
      return
      
   def get_mesh(self, path):
      """ Create or load a mesh, given a configuration specified in the simulation configuration file. 
      
      :param str path: The path to the mesh file.
      :returns: A Mesh object.
      """
      LOG.info("Creating/loading mesh...")
      
      dimension = self.options["dimension"]
      try:
         # Unit mesh whose vertex coordinates lie in the range [0, 1] along all axes.
         if(libspud.have_option("/geometry/mesh/unit_mesh")):
            number_of_nodes = libspud.get_option("/geometry/mesh/unit_mesh/number_of_nodes")
            if(dimension == 1):
               mesh = UnitIntervalMesh(number_of_nodes[0])
            elif(dimension == 2):
               mesh = UnitSquareMesh(number_of_nodes[0], number_of_nodes[1])
            elif(dimension == 3):
               mesh = UnitCubeMesh(number_of_nodes[0], number_of_nodes[1], number_of_nodes[2])
            else:
               raise ValueError("Unsupported dimension.")
               
         # Interval mesh whose vertex coordinates lie in the range [0, L].
         elif(libspud.have_option("/geometry/mesh/interval_mesh")):
            
            L = libspud.get_option("/geometry/mesh/interval_mesh/length")
            n = libspud.get_option("/geometry/mesh/interval_mesh/number_of_cells")
            mesh = IntervalMesh(n, L)
            
         # A user-defined mesh file (currently only supports Gmsh format).
         elif(libspud.have_option("/geometry/mesh/from_file")):
            absolute_path_to_config = os.path.dirname(os.path.abspath(path))
            # This is the path relative to the directory where the configuration file is stored.
            relative_path_to_mesh = libspud.get_option("/geometry/mesh/from_file/relative_path") 
            absolute_path_to_mesh = os.path.join(absolute_path_to_config, relative_path_to_mesh)
            if(not os.path.exists(absolute_path_to_mesh)):
               raise ValueError("The path to the mesh file '%s' does not exist." % absolute_path_to_mesh)
            else:
               mesh = Mesh(absolute_path_to_mesh)
            
         # Unknown mesh format.
         else:
            raise ValueError("Unsupported input mesh type.")
            
      except ValueError as e:
         LOG.exception(e)
         sys.exit()
         
      return mesh
      
   def run(self):
      """ Perform the simulation! """

      # Time-stepping parameters and constants
      LOG.info("Setting up a few constants...")
      T = self.options["T"]
      t = self.options["t"]
      theta = self.options["theta"]
      dt = self.options["dt"]
      dimension = self.options["dimension"]
      g_magnitude = self.options["g_magnitude"]

      # Get the function spaces
      U = self.function_spaces["VelocityFunctionSpace"]
      H = self.function_spaces["FreeSurfaceFunctionSpace"]
      
      # Is the Velocity field represented by a discontinous function space?
      dg = (U.ufl_element().family() == "Discontinuous Lagrange")
      
      # Weight u and h by theta to obtain the theta time-stepping scheme.
      assert(theta >= 0.0 and theta <= 1.0)
      LOG.info("Time-stepping scheme using theta = %g" % (theta))
      u_mid = (1.0 - theta) * self.u0 + theta * self.u
      h_mid = (1.0 - theta) * self.h0 + theta * self.h
         
      # The total height of the free surface.
      self.h_total = self.h_mean + self.h0

      # Non-linear approximation to the velocity
      u_nl = Function(U).assign(self.u0)
      
      # Second-order Adams-Bashforth velocity
      u_bash = (3.0/2.0)*self.u0 - (1.0/2.0)*self.u00
      
      # Simple P1 function space, to be used in the stabilisation routines (if applicable).
      P1 = FunctionSpace(self.mesh, "CG", 1)
      cellsize = CellSize(self.mesh)

      ###########################################################
      ################# Tentative velocity step #################
      ###########################################################
      
      # The collection of all the individual terms in their weak form.
      LOG.info("Constructing form...")
      F = 0

      # Mass term
      if(self.options["have_momentum_mass"]):
         LOG.debug("Momentum equation: Adding mass term...")
         M_momentum = (1.0/dt)*(inner(self.w, self.u) - inner(self.w, self.u0))*dx
         F += M_momentum
         
      # Append any Expression objects for weak BCs here.
      weak_bc_expressions = []
      
      # Advection term
      if(self.options["have_momentum_advection"]):
         LOG.debug("Momentum equation: Adding advection term...")

         if(self.options["integrate_advection_term_by_parts"]):
            outflow = (dot(self.u0, self.n) + abs(dot(self.u0, self.n)))/2.0

            A_momentum = -inner(dot(u_nl, grad(self.w)), u_bash)*dx - inner(dot(u_bash, grad(u_nl)), self.w)*dx
            A_momentum += inner(self.w, outflow*u_mid)*ds
            if(dg):
               # Only add interior facet integrals if we are dealing with a discontinous Galerkin discretisation.
               A_momentum += dot(outflow('+')*u_mid('+') - outflow('-')*u_mid('-'), jump(self.w))*dS

         else:
            A_momentum = inner(dot(grad(self.u), u_nl), self.w)*dx
         F += A_momentum

      # Viscous stress term. Note that the viscosity is kinematic (not dynamic).
      if(self.options["have_momentum_stress"]):
         LOG.debug("Momentum equation: Adding stress term...")
         
         viscosity = Function(H)
         
         # Background viscosity
         background_viscosity = Function(H).interpolate(Expression(libspud.get_option("/system/equations/momentum_equation/stress_term/scalar_field::Viscosity/value/constant")))
         viscosity.assign(background_viscosity)

         # Eddy viscosity
         if(self.options["have_turbulence_parameterisation"]):
            LOG.debug("Momentum equation: Adding turbulence parameterisation...")
            base_option_path = "/system/equations/momentum_equation/turbulence_parameterisation"
            # Large eddy simulation (LES)
            if(libspud.have_option(base_option_path + "/les")):
               les = LES(self.mesh, H)
               density = Constant(1.0) # We divide through by density in the momentum equation, so just set this to 1.0 for now.
               smagorinsky_coefficient = Constant(libspud.get_option(base_option_path + "/les/smagorinsky/smagorinsky_coefficient"))
               
               eddy_viscosity = Function(H)
               eddy_viscosity_lhs, eddy_viscosity_rhs = les.eddy_viscosity(u_mid, density, smagorinsky_coefficient)
               eddy_viscosity_problem = LinearVariationalProblem(eddy_viscosity_lhs, eddy_viscosity_rhs, eddy_viscosity, bcs=[])
               eddy_viscosity_solver = LinearVariationalSolver(eddy_viscosity_problem)
               
            # Add on eddy viscosity
            viscosity += eddy_viscosity

         # Stress tensor: tau = grad(u) + transpose(grad(u)) - (2/3)*div(u)
         if(not dg):
            # Perform a double dot product of the stress tensor and grad(w).
            K_momentum = -viscosity*inner(grad(self.u) + grad(self.u).T, grad(self.w))*dx
            K_momentum += viscosity*(2.0/3.0)*inner(div(self.u)*Identity(dimension), grad(self.w))*dx
         else:
            # Interior penalty method
            cellsize = Constant(0.2) # In general, we should use CellSize(self.mesh) instead.
            alpha = 1/cellsize # Penalty parameter.
            
            K_momentum = -viscosity('+')*inner(grad(u_mid), grad(self.w))*dx
            for dim in range(self.options["dimension"]):
               K_momentum += -viscosity('+')*(alpha('+')/cellsize('+'))*dot(jump(self.w[dim], self.n), jump(u_mid[dim], self.n))*dS
               K_momentum += viscosity('+')*dot(avg(grad(self.w[dim])), jump(u_mid[dim], self.n))*dS + viscosity('+')*dot(jump(self.w[dim], self.n), avg(grad(u_mid[dim])))*dS

         F -= K_momentum # Negative sign here because we are bringing the stress term over from the RHS.

      # The gradient of the height of the free surface, h
      LOG.debug("Momentum equation: Adding gradient term...")
      C_momentum = -g_magnitude*inner(self.w, grad(self.h0))*dx
      F -= C_momentum

      # Quadratic drag term in the momentum equation
      if(self.options["have_drag"]):
         LOG.debug("Momentum equation: Adding drag term...")
         
         base_option_path = "/system/equations/momentum_equation/drag_term"
         
         # Get the bottom drag/friction coefficient.
         LOG.debug("Momentum equation: Adding bottom drag contribution...")
         drag_expression = ExpressionFromOptions(path=base_option_path+"/scalar_field::BottomDragCoefficient/value", t=t).get_expression()
         drag_coefficient = Function(H).interpolate(drag_expression)
         
         # Add on the turbine drag, if provided.
         self.array = None
         
         # Magnitude of the velocity field
         magnitude = sqrt(dot(self.u0, self.u0))
         
         # Form the drag term
         D_momentum = -inner(self.w, (drag_coefficient*magnitude/self.h_total)*u_nl)*dx
         F -= D_momentum

      # Add in any source terms
      if(self.options["have_momentum_source"]):
         LOG.debug("Momentum equation: Adding source term...")
         momentum_source_expression = ExpressionFromOptions(path = "/system/equations/momentum_equation/source_term/vector_field::Source/value", t=t).get_expression()
         momentum_source_function = Function(U)
         F -= inner(self.w, momentum_source_function.interpolate(momentum_source_expression))*dx
         
      # Add in any SU stabilisation
      if(self.options["have_su_stabilisation"]):
         LOG.debug("Momentum equation: Adding streamline-upwind stabilisation term...")
         stabilisation = Stabilisation(self.mesh, P1, cellsize)
         
         magnitude = magnitude_vector(self.u0, P1)

         # Bound the values for the magnitude below by 1.0e-9 for numerical stability reasons.
         u_nodes = magnitude.vector()
         near_zero = numpy.array([1.0e-9 for i in range(len(u_nodes))])
         u_nodes.set_local(numpy.maximum(u_nodes.array(), near_zero))

         diffusivity = ExpressionFromOptions(path = "/system/equations/momentum_equation/stress_term/scalar_field::Viscosity/value", t=self.options["t"]).get_expression()
         diffusivity = Function(H).interpolate(diffusivity) # Background viscosity
         grid_pe = grid_peclet_number(diffusivity, magnitude, P1, cellsize)
   
         # Bound the values for grid_pe below by 1.0e-9 for numerical stability reasons. 
         grid_pe_nodes = grid_pe.vector()
         values = numpy.array([1.0e-9 for i in range(len(grid_pe_nodes))])
         grid_pe_nodes.set_local(numpy.maximum(grid_pe_nodes.array(), values))

         F += stabilisation.streamline_upwind(self.w, self.u0, magnitude, grid_pe)

      LOG.info("Form construction complete.")

      ##########################################################
      ################ Pressure correction step ################
      ##########################################################   
      u_tent = Function(U)
      u_tent_nl = theta*u_tent + (1.0-theta)*self.u0
      F_h_corr = inner(self.v, (self.h - self.h0))*dx \
                 + g_magnitude*(dt**2)*(theta**2)*self.h_total*inner(grad(self.v), grad(self.h - self.h0))*dx \
                 + dt*self.v*div(self.h_total*u_tent_nl)*dx
                 
      ##########################################################
      ################ Velocity correction step ################
      ##########################################################   
      h1 = Function(H)
      u1 = Function(U)
      F_u_corr = (1.0/dt)*inner(self.w, self.u - u_tent)*dx + g_magnitude*theta*inner(self.w, grad(h1 - self.h0))*dx
                 
      LOG.info("Applying strong Dirichlet boundary conditions...")
      # Get all the Dirichlet boundary conditions for the Velocity field
      bcs_u = []; bcs_u2 = []
      bcs_h = []
      bc_expressions = []
      for i in range(0, libspud.option_count("/system/core_fields/vector_field::Velocity/boundary_condition")):
         if(libspud.have_option("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet" % i) and
            not libspud.have_option("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet/apply_weakly" % i)):
            expr = ExpressionFromOptions(path = ("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet" % i), t=t).get_expression()
            # Surface IDs on the domain boundary
            surface_ids = libspud.get_option("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/surface_ids" % i)
            method = ("geometric" if dg else "topological")
            bc = DirichletBC(U, expr, surface_ids, method=method)
            bcs_u.append(bc)
            bc_expressions.append(expr)
            LOG.debug("Applying Velocity BC #%d strongly to surface IDs: %s" % (i, surface_ids))

      for i in range(0, libspud.option_count("/system/core_fields/vector_field::Velocity/boundary_condition")):
         if(libspud.have_option("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet" % i) and
            not libspud.have_option("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet/apply_weakly" % i)):
            expr = ExpressionFromOptions(path = ("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/type::dirichlet" % i), t=t).get_expression()
            # Surface IDs on the domain boundary
            surface_ids = libspud.get_option("/system/core_fields/vector_field::Velocity/boundary_condition[%d]/surface_ids" % i)
            method = ("geometric" if dg else "topological")
            bc = DirichletBC(U, expr, surface_ids, method=method)
            bcs_u2.append(bc)
            bc_expressions.append(expr)
            LOG.debug("Applying Velocity BC #%d strongly to surface IDs: %s" % (i, surface_ids))


      for i in range(0, libspud.option_count("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition/type::dirichlet")):
         if(libspud.have_option("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/type::dirichlet" % i) and
            not(libspud.have_option("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/type::dirichlet/apply_weakly" % i))):
            expr = ExpressionFromOptions(path = ("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/type::dirichlet" % i), t=t).get_expression()
            # Surface IDs on the domain boundary
            surface_ids = libspud.get_option("/system/core_fields/scalar_field::FreeSurfacePerturbation/boundary_condition[%d]/surface_ids" % i)
            method = ("geometric" if dg else "topological")
            bc = DirichletBC(H, expr, surface_ids, method=method)
            bcs_h.append(bc)
            bc_expressions.append(expr)
            LOG.debug("Applying FreeSurfacePerturbation BC #%d strongly to surface IDs: %s" % (i, surface_ids))
            
      # Prepare solver_parameters dictionary
      LOG.debug("Defining solver_parameters dictionary...")
      solver_parameters = {'ksp_monitor': True, 'ksp_view': False, 'pc_view': False, 'snes_type': 'ksponly', 'ksp_max_it':10000} # NOTE: use 'snes_type': 'newtonls' for production runs.
      
      # KSP (iterative solver) options
      solver_parameters["ksp_type"] = libspud.get_option("/system/solver/iterative_method/name")
      solver_parameters["ksp_rtol"] = libspud.get_option("/system/solver/relative_error")
      
      solver_parameters['ksp_converged_reason'] = True
      solver_parameters['ksp_monitor_true_residual'] = True
      
      # Preconditioner options
      solver_parameters["pc_type"] = libspud.get_option("/system/solver/preconditioner/name")
      # Fieldsplit sub-options
      if(solver_parameters["pc_type"] == "fieldsplit"):
         LOG.debug("Setting up the fieldsplit preconditioner...")
         solver_parameters["pc_fieldsplit_type"] = libspud.get_option("/system/solver/preconditioner::fieldsplit/type/name")
         if(solver_parameters["pc_fieldsplit_type"] == "schur"):
            solver_parameters["pc_fieldsplit_schur_fact_type"] = libspud.get_option("/system/solver/preconditioner::fieldsplit/type::schur/fact_type/name")
         solver_parameters["fieldsplit_0_ksp_type"] = libspud.get_option("/system/solver/preconditioner::fieldsplit/block_0_ksp_type/iterative_method/name")
         solver_parameters["fieldsplit_1_ksp_type"] = libspud.get_option("/system/solver/preconditioner::fieldsplit/block_1_ksp_type/iterative_method/name")

         if(libspud.get_option("/system/solver/preconditioner::fieldsplit/block_0_pc_type/preconditioner/name") != "ilu"):
            solver_parameters["fieldsplit_0_pc_type"] = libspud.get_option("/system/solver/preconditioner::fieldsplit/block_0_pc_type/preconditioner/name")
            solver_parameters["fieldsplit_1_pc_type"] = libspud.get_option("/system/solver/preconditioner::fieldsplit/block_1_pc_type/preconditioner/name")
            
         # Enable inner iteration monitors.
         solver_parameters["fieldsplit_0_ksp_monitor"] = True
         solver_parameters["fieldsplit_1_ksp_monitor"] = True
         solver_parameters["fieldsplit_0_pc_factor_shift_type"] = 'INBLOCKS'
         solver_parameters["fieldsplit_1_pc_factor_shift_type"] = 'INBLOCKS'
         
      # Construct the solver objects
      problem_tent = LinearVariationalProblem(lhs(F), rhs(F), u_tent, bcs=bcs_u)
      solver_tent = LinearVariationalSolver(problem_tent, solver_parameters={'ksp_monitor': True, 
                                                                              'ksp_view': False, 
                                                                              'pc_view': False, 
                                                                              'pc_type': 'sor',
                                                                              'ksp_type': 'gmres',
                                                                              'ksp_rtol': 1.0e-7})

      problem_h_corr = LinearVariationalProblem(lhs(F_h_corr), rhs(F_h_corr), h1, bcs=bcs_h)
      solver_h_corr = LinearVariationalSolver(problem_h_corr, solver_parameters={'ksp_monitor': True, 
                                                                                 'ksp_view': False, 
                                                                                 'pc_view': False, 
                                                                                 'pc_type': 'sor',
                                                                                 'ksp_type': 'gmres',
                                                                                 'ksp_rtol': 1.0e-7})
      
     
      problem_u_corr = LinearVariationalProblem(lhs(F_u_corr), rhs(F_u_corr), u1, bcs=bcs_u2)
      solver_u_corr = LinearVariationalSolver(problem_u_corr, solver_parameters={'ksp_monitor': True, 
                                                                                 'ksp_view': False, 
                                                                                 'pc_view': False, 
                                                                                 'pc_type': 'sor',
                                                                                 'ksp_type': 'gmres',
                                                                                 'ksp_rtol': 1.0e-7})
      
      
      t += dt
      iterations_since_dump = 1
      iterations_since_checkpoint = 1
      
      # PETSc solver run-times
      from petsc4py import PETSc
      main_solver_stage = PETSc.Log.Stage('Main block-coupled system solve')

      total_solver_time = 0.0
      # The time-stepping loop
      LOG.info("Entering the time-stepping loop...")
      EPSILON = 1.0e-14
      
      while t <= T + EPSILON: # A small value EPSILON is added here in case of round-off error.
         LOG.info("t = %g" % t)

         while True:
            ## Update any time-dependent Functions and Expressions.
            
            # Re-compute the velocity magnitude and grid Peclet number fields.
            if(self.options["have_su_stabilisation"]):
               magnitude.assign(magnitude_vector(self.u0, P1))

               # Bound the values for the magnitude below by 1.0e-9 for numerical stability reasons.
               u_nodes = magnitude.vector()
               near_zero = numpy.array([1.0e-9 for i in range(len(u_nodes))])
               u_nodes.set_local(numpy.maximum(u_nodes.array(), near_zero))
               
               grid_pe.assign(grid_peclet_number(diffusivity, magnitude, P1, cellsize))
         
               # Bound the values for grid_pe below by 1.0e-9 for numerical stability reasons. 
               grid_pe_nodes = grid_pe.vector()
               values = numpy.array([1.0e-9 for i in range(len(grid_pe_nodes))])
               grid_pe_nodes.set_local(numpy.maximum(grid_pe_nodes.array(), values))

            if(self.options["have_turbulence_parameterisation"]):
               eddy_viscosity_solver.solve()
               viscosity.assign(background_viscosity + eddy_viscosity)

            # Time-dependent source terms
            if(self.options["have_momentum_source"]):
               momentum_source_expression.t = t
               momentum_source_function.interpolate(momentum_source_expression)
            if(self.options["have_continuity_source"]):
               continuity_source_expression.t = t
               continuity_source_function.interpolate(continuity_source_expression)
               
            # Update any time-varying DirichletBC objects.
            for expr in bc_expressions:
               expr.t = t
            for expr in weak_bc_expressions:
               expr.t = t
                      
            # Solve the system of equations!
            start_solver_time = mpi4py.MPI.Wtime()
            main_solver_stage.push()
            LOG.debug("Solving the system of equations...")
            
            solver_tent.solve()
            solver_h_corr.solve()
            solver_u_corr.solve()
            
            main_solver_stage.pop()
            end_solver_time = mpi4py.MPI.Wtime()
            total_solver_time += (end_solver_time - start_solver_time)

            # Move to next time step  
            if(steady_state(u1, u_nl, 1e-7)):
               break
            u_nl.assign(u1)
      

         self.u00.assign(self.u0)
         self.u0.assign(u1)
         self.h0.assign(h1)
         t += dt
         iterations_since_dump += 1
         iterations_since_checkpoint += 1
         LOG.debug("Moving to next time level...")

         # Write the solution to file.
         if((self.options["dump_period"] is not None) and (dt*iterations_since_dump >= self.options["dump_period"])):
            LOG.debug("Writing data to file...")
            self.output_functions["Velocity"].assign(u1)
            self.output_files["Velocity"] << self.output_functions["Velocity"]
            self.output_functions["FreeSurfacePerturbation"].assign(h1)
            self.output_files["FreeSurfacePerturbation"] << self.output_functions["FreeSurfacePerturbation"]
            iterations_since_dump = 0 # Reset the counter.

         # Print out the total power generated by turbines.
         if(self.options["have_drag"] and self.array is not None):
            LOG.info("Power = %.2f" % self.array.power(u1, density=1000))
            
         # Checkpointing
         if((self.options["checkpoint_period"] is not None) and (dt*iterations_since_checkpoint >= self.options["checkpoint_period"])):
            LOG.debug("Writing checkpoint data to file...")
            self.solution.dat.save("checkpoint")
            iterations_since_checkpoint = 0 # Reset the counter.
            
         # Check whether a steady-state has been reached.
         if(steady_state(u1, self.u0, self.options["steady_state_tolerance"]) and steady_state(h1, self.h0, self.options["steady_state_tolerance"])):
            LOG.info("Steady-state attained. Exiting the time-stepping loop...")
            break

         self.compute_diagnostics()            


      LOG.info("Out of the time-stepping loop.")

      LOG.debug("Total solver time: %.2f" % (total_solver_time))

      return

if(__name__ == "__main__"):
   import signal
   import argparse
   
   # Get the Git revision of the code.
   revision = get_git_revision(os.path.dirname(os.path.realpath(__file__)))
   if(revision is not None):
      LOG.info("Git revision: %s" % revision)

   # Parse options and arguments from the command line
   LOG.info("Parsing command line arguments...")
   usage = "Usage: python /path/to/shallow_water.py [options] path/to/simulation_setup_file.swml"
   parser = argparse.ArgumentParser(description="The shallow water model in the Firedrake-Fluids CFD code.")
   parser.add_argument("-c", "--checkpoint", action="store", default=None, type=str, help="Initialise field values from a specified checkpoint file.", metavar="CHECKPOINT_FILE")
   parser.add_argument("path", help="The path to the simulation configuration file (with a .swml extension).", action="store", type=str)
   args = parser.parse_args()
        
   # Exit if SIGINT is detected.
   signal.signal(signal.SIGINT, signal.SIG_DFL)
   
   if(os.path.exists(args.path)):
   
      simulation_start_time = mpi4py.MPI.Wtime()
      
      # Set up a shallow water simulation.
      sw = ShallowWater(path=args.path, checkpoint=args.checkpoint)
      
      # Solve the shallow water equations!
      sw.run()
      
      simulation_end_time = mpi4py.MPI.Wtime()
      
      LOG.info("Total simulation run-time = %.2f s" % (simulation_end_time - simulation_start_time))
      
   else:
      LOG.error("The path to the simulation setup file does not exist.")
      sys.exit()

