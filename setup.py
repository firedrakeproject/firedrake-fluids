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

from distutils.core import setup

setup(name='Firedrake-Fluids',
      version='0.2-dev',
      description='A collection of numerical models for fluid flow simulation, using the Firedrake framework for the portable solution of the underlying model equations.',
      author='Imperial College London',
      url='https://github.com/firedrakeproject/firedrake-fluids',
      packages=['firedrake_fluids'],
      package_dir = {'firedrake_fluids': 'firedrake_fluids'}
     )

