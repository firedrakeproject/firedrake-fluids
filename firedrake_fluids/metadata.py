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
