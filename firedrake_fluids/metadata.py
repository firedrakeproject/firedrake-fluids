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

import subprocess
from firedrake_fluids import LOG

def get_git_revision(cwd=None):
   """ Return the Git HEAD revision in the form of a SHA-1 hash. """
   # Adapted from the pybench code by Florian Rathgeber.
   # https://github.com/firedrakeproject/pybench/blob/master/pybench.py
   try:
      revision = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=cwd).strip()
   except OSError as e:
      LOG.exception(e)
      revision = None
   return revision
