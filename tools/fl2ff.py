#!/usr/bin/env python

#    Copyright (C) 2014 Imperial College London.

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

import libspud

class VelocityBoundaryCondition:
   def __init__(self, option_path):
      print option_path
      self.name = libspud.get_option(option_path + "/name")
      self.surface_ids = libspud.get_option(option_path + "/surface_ids")
      self.type = libspud.get_option(option_path + "/type/name")
      
      if(self.type == "dirichlet"):
         self.weak = libspud.have_option(option_path + "/type::%s/apply_weakly" % self.type)
         self.value = []
         R = ["x", "y", "z"]
         for r in R:
            if(libspud.have_option(option_path + "/type::%s/align_bc_with_cartesian/%s_component" % (self.type, r))):
               c = libspud.get_child_name(option_path + "/type::%s/align_bc_with_cartesian/%s_component/" % (self.type), 1)
               self.value.append(libspud.get_option(option_path + "/type::%s/align_bc_with_cartesian/%s_component/%s" % (self.type, r, c)))
            else:
               self.value.append(None)
         print self.value
      else:
         self.weak = None
         self.value = [None, None, None]

class PressureBoundaryCondition:
   def __init__(self, option_path):
      print option_path
      self.name = libspud.get_option(option_path + "/name")
      self.surface_ids = libspud.get_option(option_path + "/surface_ids")
      self.type = libspud.get_option(option_path + "/type/name")
      
      if(self.type == "dirichlet"):
         c = libspud.get_child_name(option_path + "/type::%s/", 1)
         self.value.append(libspud.get_option(option_path + "/type::%s/align_bc_with_cartesian/%s_component/%s" % (self.type, c)))
      else:
         self.value = None
         
class FunctionSpace:
   def __init__(self, option_path, index):
      self.name = libspud.get_option(option_path + "/mesh[%d]/name" % index)
      self.degree = libspud.get_option(option_path + "/mesh[%d]/from_mesh/mesh_shape/polynomial_degree" % index)
      if(libspud.have_option(option_path + "/mesh[%d]/from_mesh/mesh_continuity" % index)):
         if(libspud.get_option(option_path + "/mesh[%d]/from_mesh/mesh_continuity" % index) == "continuous"):
            self.family = "Continuous Lagrange"
         else:
            self.family = "Discontinuous Lagrange"
      else:
         self.family = "Continuous Lagrange"
         
def convert(path):
   
   # Read in Fluidity simulation options
   
   libspud.clear_options()
   libspud.load_options(path)
   
   # Simulation name
   simulation_name = libspud.get_option("/simulation_name")
   
   # Geometry
   base = "/geometry"
   dimension = libspud.get_option(base + "/dimension")
   mesh_path = libspud.get_option(base + "/mesh::CoordinateMesh/from_file/file_name") + ".msh" # FIXME: Always assumes gmsh format.
  
   # Function spaces
   velocity_function_space = FunctionSpace("/geometry", 1)
   freesurface_function_space = FunctionSpace("/geometry", 2)
   
   # I/O
   base = "/io"
   dump_format = libspud.get_option(base + "/dump_format")
   dump_period = libspud.get_option(base + "/dump_period/constant")
      
   # Timestepping
   base = "/timestepping"
   current_time = libspud.get_option(base + "/current_time")
   timestep = libspud.get_option(base + "/timestep")
   finish_time = libspud.get_option(base + "/finish_time")

   ## Steady-state
   if(libspud.have_option(base + "/steady_state")):
      if(libspud.have_option(base + "/steady_state/tolerance")):
         steady_state = libspud.get_option(base + "/steady_state/tolerance")
      else:
         steady_state = 1e-7
   else:
      steady_state = None   
   
   # Gravity
   g_magnitude = libspud.get_option("/physical_parameters/gravity/magnitude")
   
   # Velocity field (momentum equation)
   base = "/material_phase[0]/vector_field::Velocity"
   
   ## Depth (free surface mean height)
   c = libspud.get_child_name(base + "/prognostic/equation::ShallowWater/scalar_field::BottomDepth/prescribed/value::WholeMesh/", 1)
   depth = libspud.get_option(base + "/prognostic/equation::ShallowWater/scalar_field::BottomDepth/prescribed/value::WholeMesh/%s" % c)
   
   ## Bottom drag coefficient
   if(libspud.have_option(base + "/prognostic/equation::ShallowWater/bottom_drag")):
      c = libspud.get_child_name(base + "/prognostic/equation::ShallowWater/bottom_drag/scalar_field::BottomDragCoefficient/prescribed/value::WholeMesh/", 1)
      bottom_drag = libspud.get_option(base + "/prognostic/equation::ShallowWater/bottom_drag/scalar_field::BottomDragCoefficient/prescribed/value::WholeMesh/%s" % c)
   else:
      bottom_drag = None
      
   ## Viscosity
   if(libspud.have_option(base + "/prognostic/tensor_field::Viscosity")):
      viscosity = libspud.get_option(base + "/prognostic/tensor_field::Viscosity/prescribed/value::WholeMesh/anisotropic_symmetric/constant")[0][0]
   else:
      viscosity = None

   ## Momentum source
   if(libspud.have_option(base + "/prognostic/vector_field::Source")):
      c = libspud.get_child_name(base + "/prognostic/vector_field::Source/prescribed/value::WholeMesh/", 1)
      momentum_source = libspud.get_option(base + "/prognostic/vector_field::Source/prescribed/value::WholeMesh/%s" % c)
   else:
      momentum_source = None
      
   ## Initial condition
   if(libspud.have_option(base + "/prognostic/initial_condition::WholeMesh")):
      c = libspud.get_child_name(base + "/prognostic/initial_condition::WholeMesh/", 1)
      velocity_initial_condition = libspud.get_option(base + "/prognostic/initial_condition::WholeMesh/%s" % c)
   else:
      velocity_initial_condition = 0.0
      
   ## Boundary conditions
   number_of_bcs = libspud.option_count(base + "/prognostic/boundary_conditions")
   velocity_bcs = []
   for i in range(number_of_bcs):
      velocity_bcs.append(VelocityBoundaryCondition(base + "/prognostic/boundary_conditions[%d]" % i))
   
   
   # Pressure field (continuity equation)
   base = "/material_phase[0]/scalar_field::Pressure"
   integrate_by_parts = libspud.have_option(base + "/prognostic/spatial_discretisation/continuous_galerkin/integrate_continuity_by_parts")
   
   ## Initial condition
   if(libspud.have_option(base + "/prognostic/initial_condition::WholeMesh")):
      c = libspud.get_child_name(base + "/prognostic/initial_condition::WholeMesh/", 1)
      pressure_initial_condition = libspud.get_option(base + "/prognostic/initial_condition::WholeMesh/%s" % c)
   else:
      pressure_initial_condition = 0.0
      
   ## Boundary conditions
   number_of_bcs = libspud.option_count(base + "/prognostic/boundary_conditions")
   pressure_bcs = []
   for i in range(number_of_bcs):
      pressure_bcs.append(PressureBoundaryCondition(base + "/prognostic/boundary_conditions[%d]" % i))
   
   ## Continuity source
   if(libspud.have_option(base + "/prognostic/scalar_field::Source")):
      c = libspud.get_child_name(base + "/prognostic/scalar_field::Source/prescribed/value::WholeMesh/", 1)
      continuity_source = libspud.get_option(base + "/prognostic/scalar_field::Source/prescribed/value::WholeMesh/%s" % c)
   else:
      continuity_source = None
   





   
   
   # Write out to a Firedrake-Fluids simulation configuration file
   libspud.clear_options()
   
   # Create a bare-bones .swml file to add to.
   f = open("dummy.swml", "w")
   f.write("<?xml version='1.0' encoding='utf-8'?>\n")
   f.write("<shallow_water_options>\n")
   f.write("</shallow_water_options>\n")
   f.close()
   
   libspud.load_options("dummy.swml")

   # Simulation name
   libspud.set_option("/simulation_name", simulation_name)
   
   # Geometry
   base = "/geometry"
   libspud.set_option(base + "/dimension", dimension)
   libspud.set_option(base + "/mesh/from_file/relative_path", mesh_path)
   
   # Function spaces
   base = "/function_spaces"
   libspud.set_option(base + "/function_space::VelocityFunctionSpace/degree", velocity_function_space.degree)
   libspud.set_option(base + "/function_space::VelocityFunctionSpace/family", velocity_function_space.family)
   libspud.set_option(base + "/function_space::FreeSurfaceFunctionSpace/degree", freesurface_function_space.degree)
   libspud.set_option(base + "/function_space::FreeSurfaceFunctionSpace/family", freesurface_function_space.family)
   
   # I/O
   base = "/io"
   libspud.set_option(base + "/dump_format", dump_format)
   libspud.set_option(base + "/dump_period", dump_period)
   
   # Timestepping
   base = "/timestepping"
   print timestep
   libspud.set_option(base + "/current_time", current_time)
   try:
      libspud.set_option(base + "/timestep", timestep)
   except:
      pass
   libspud.set_option(base + "/finish_time", finish_time)
   
   ## Steady-state
   if(steady_state):
      libspud.set_option(base + "/steady_state/tolerance", steady_state)
      
   # Gravity
   libspud.set_option("/physical_parameters/gravity/magnitude", g_magnitude)
   
   # System/Core Fields: Velocity
   base = "/system/core_fields/vector_field::Velocity"
   
   ## Initial condition
   if(isinstance(velocity_initial_condition, str)):
      libspud.set_option(base + "/initial_condition/python", velocity_initial_condition)
   else:
      libspud.set_option(base + "/initial_condition/constant", velocity_initial_condition)
   
   ## Boundary conditions
   try:
      for i in range(len(velocity_bcs)):
         libspud.set_option(base + "/boundary_condition::%s/surface_ids" % velocity_bcs[i].name, velocity_bcs[i].surface_ids)
         libspud.set_option_attribute(base + "/boundary_condition::%s/type/name" % velocity_bcs[i].name, velocity_bcs[i].type)
         
         if(velocity_bcs[i].type == "dirichlet"):
            if(isinstance(velocity_bcs[i].value, str)):
               libspud.set_option(base + "/boundary_condition::%s/type::dirichlet/value/python" % velocity_bcs[i].name, velocity_bcs[i].value)
            else:
               libspud.set_option(base + "/boundary_condition::%s/type::dirichlet/value/constant" % velocity_bcs[i].name, velocity_bcs[i].value)
   except:
      pass

   # System/Core Fields: FreeSurfacePerturbation
   base = "/system/core_fields/scalar_field::FreeSurfacePerturbation"
   
   #FIXME: Pressure initial and boundary conditions are multiplied by 'g' in Fluidity, but not in Firedrake-Fluids.
   ## Initial condition
   if(isinstance(pressure_initial_condition, str)):
      libspud.set_option(base + "/initial_condition/python", pressure_initial_condition)
   else:
      libspud.set_option(base + "/initial_condition/constant", pressure_initial_condition)
   
   ## Boundary conditions
   try:
      for i in range(len(pressure_bcs)):
         libspud.set_option(base + "/boundary_condition::%s/surface_ids" % pressure_bcs[i].name, pressure_bcs[i].surface_ids)
         libspud.set_option(base + "/boundary_condition::%s/type/name" % pressure_bcs[i].name, pressure_bcs[i].type)
         
         if(pressure_bcs[i].type == "dirichlet"):
            if(isinstance(pressure_bcs[i].value, str)):
               libspud.set_option(base + "/boundary_condition::%s/type::dirichlet/value/python" % pressure_bcs[i].name, pressure_bcs[i].value)
            else:
               libspud.set_option(base + "/boundary_condition::%s/type::dirichlet/value/constant" % pressure_bcs[i].name, pressure_bcs[i].value)
   except:
      pass

   # System/Core Fields: FreeSurfaceMean
   base = "/system/core_fields/scalar_field::FreeSurfaceMean"
   if(isinstance(depth, str)):
      libspud.set_option(base + "/value/python", depth)
   else:
      libspud.set_option(base + "/value/constant", depth)
      
   
   # Equations: Continuity equation
   base = "/system/equations/continuity_equation"
   libspud.set_option(base + "/integrate_by_parts", integrate_by_parts)

   ## Source term
   if(continuity_source is not None):
      if(isinstance(continuity_source, str)):
         libspud.set_option(base + "/source_term/scalar_field::Source/value/python", continuity_source)
      else:
         libspud.set_option(base + "/source_term/scalar_field::Source/value/constant", continuity_source)
         
   # Equations: Momentum equation
   base = "/system/equations/momentum_equation"
   
   ## Viscosity
   if(viscosity is not None):
      if(isinstance(viscosity, str)):
         libspud.set_option(base + "/stress_term/scalar_field::Viscosity/value/python", viscosity)
      else:
         libspud.set_option(base + "/stress_term/scalar_field::Viscosity/value/constant", viscosity)

   ## Bottom drag
   if(bottom_drag is not None):
      if(isinstance(bottom_drag, str)):
         libspud.set_option(base + "/drag_term/scalar_field::BottomDragCoefficient/value/python", bottom_drag)
      else:
         libspud.set_option(base + "/drag_term/scalar_field::BottomDragCoefficient/value/constant", bottom_drag)
   
   ## Source term
   if(momentum_source is not None):
      if(isinstance(momentum_source, str)):
         libspud.set_option(base + "/source_term/vector_field::Source/value/python", momentum_source)
      else:
         libspud.set_option(base + "/source_term/vector_field::Source/value/constant", momentum_source)
   
   
   # Write all the applied options to file.
   libspud.write_options("dummy.swml")
   
   return

if(__name__ == "__main__"):
   import argparse

   # Parse options and arguments from the command line
   usage = "Usage: python fl2ff.py [options] path/to/fluidity_setup_file.flml path/to/firedrake-fluids_setup_file.swml"
   parser = argparse.ArgumentParser(description="Converts a Fluidity shallow water setup file to a Firedrake-Fluids shallow water setup file.")
   parser.add_argument("path_to_fluidity_setup_file", help="The path to the Fluidity simulation configuration file (with a .flml extension).", action="store", type=str)
   parser.add_argument("path_to_firedrake-fluids_setup_file", help="The path to the Firedrake-Fluids simulation configuration file (with a .swml extension).", action="store", type=str)
   args = parser.parse_args()

   convert(args.path_to_fluidity_setup_file)

