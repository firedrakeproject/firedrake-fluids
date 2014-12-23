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
               if(libspud.have_option(option_path + "/type::%s/align_bc_with_cartesian/%s_component/constant" % (self.type, r))):
                  self.value.append(libspud.get_option(option_path + "/type::%s/align_bc_with_cartesian/%s_component/constant" % (self.type, r)))
               else:
                  self.value.append(libspud.get_option(option_path + "/type::%s/align_bc_with_cartesian/%s_component/python" % (self.type, r)))
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
         if(libspud.have_option(option_path + "/type::%s/constant" % self.type)):
            self.value.append(libspud.get_option(option_path + "/type::%s/align_bc_with_cartesian/%s_component/constant" % self.type))
         else:
            self.value.append(libspud.get_option(option_path + "/type::%s/align_bc_with_cartesian/%s_component/python" % self.type))
      else:
         self.value = None
         
class FunctionSpace:
   def __init__(self, option_path, index):
      self.name = libspud.get_option(option_path + "/mesh[%d]/name" % index)

      self.degree = libspud.get_option(option_path + "/mesh[%d]/from_mesh/mesh_shape/polynomial_degree" % index)
      
      if(libspud.have_option(option_path + "/mesh[%d]/from_mesh/mesh_continuity" % index)):
         if(libspud.get_option(option_path + "/mesh[%d]/from_mesh/mesh_continuity" % index) == "continuous"):
            self.continuity = "Continuous Lagrangian"
         else:
            self.continuity = "Discontinuous Lagrangian"
      else:
         self.continuity = "Continuous Lagrangian"
         
def convert(path):
   
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
   if(libspud.have_option(base + "/prognostic/equation::ShallowWater/scalar_field::BottomDepth/prescribed/value::WholeMesh/constant")):
      depth = libspud.get_option(base + "/prognostic/equation::ShallowWater/scalar_field::BottomDepth/prescribed/value::WholeMesh/constant")
   else:
      depth = libspud.get_option(base + "/prognostic/equation::ShallowWater/scalar_field::BottomDepth/prescribed/value::WholeMesh/python")
   
   ## Bottom drag coefficient
   if(libspud.have_option(base + "/prognostic/equation::ShallowWater/bottom_drag")):
      if(libspud.have_option(base + "/prognostic/equation::ShallowWater/bottom_drag/scalar_field::BottomDragCoefficient/prescribed/value::WholeMesh/constant")):
         bottom_drag = libspud.get_option(base + "/prognostic/equation::ShallowWater/bottom_drag/scalar_field::BottomDragCoefficient/prescribed/value::WholeMesh/constant")
      else:
         bottom_drag = libspud.get_option(base + "/prognostic/equation::ShallowWater/bottom_drag/scalar_field::BottomDragCoefficient/prescribed/value::WholeMesh/python")
   else:
      bottom_drag = None
      
   ## Viscosity
   if(libspud.have_option(base + "/prognostic/tensor_field::Viscosity")):
      viscosity = libspud.get_option(base + "/prognostic/tensor_field::Viscosity/prescribed/value::WholeMesh/anisotropic_symmetric/constant")[0][0]
   else:
      viscosity = None
   
   
   ## Initial condition
   if(libspud.have_option(base + "/prognostic/initial_condition::WholeMesh")):
      if(libspud.have_option(base + "/prognostic/initial_condition::WholeMesh/constant")):
         initial_condition = libspud.get_option(base + "/prognostic/initial_condition::WholeMesh/constant")
      else:
         initial_condition = libspud.get_option(base + "/prognostic/initial_condition::WholeMesh/python")
   else:
      initial_condition = 0.0
      
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
      if(libspud.have_option(base + "/prognostic/initial_condition::WholeMesh/constant")):
         initial_condition = libspud.get_option(base + "/prognostic/initial_condition::WholeMesh/constant")
      else:
         initial_condition = libspud.get_option(base + "/prognostic/initial_condition::WholeMesh/python")
   else:
      initial_condition = 0.0
      
   ## Boundary conditions
   number_of_bcs = libspud.option_count(base + "/prognostic/boundary_conditions")
   pressure_bcs = []
   for i in range(number_of_bcs):
      pressure_bcs.append(PressureBoundaryCondition(base + "/prognostic/boundary_conditions[%d]" % i))
   
   
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

