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
import logging

def get_fluidity_options(path_to_fluidity_setup_file):
   
   return

if(__name__ == "__main__"):
   import argparse

   # Parse options and arguments from the command line
   logging.info("Parsing command line arguments...")
   usage = "Usage: python fl2ff.py [options] path/to/fluidity_setup_file.flml path/to/firedrake-fluids_setup_file.swml"
   parser = argparse.ArgumentParser(description="Converts a Fluidity shallow water setup file to a Firedrake-Fluids shallow water setup file.")
   parser.add_argument("path_to_fluidity_setup_file", help="The path to the Fluidity simulation configuration file (with a .flml extension).", action="store", type=str)
   parser.add_argument("path_to_firedrake-fluids_setup_file", help="The path to the Firedrake-Fluids simulation configuration file (with a .swml extension).", action="store", type=str)
   args = parser.parse_args()

