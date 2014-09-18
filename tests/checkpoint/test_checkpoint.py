""" Tests the checkpointing functionality. """

import os
import pytest
import numpy
from firedrake import *

from firedrake_fluids.shallow_water import *

cwd = os.path.dirname(os.path.abspath(__file__))

def checkpoint():

   # Initialise the simulation from scratch.
   sw = ShallowWater(path=os.path.join(cwd, "checkpoint.swml"), checkpoint=None)
   ux_initial = sw.solution_old.split()[0]
   h_initial = sw.solution_old.split()[1]

   # Initialise the simulation using the data in checkpoint.npy for the initial conditions.
   sw = ShallowWater(path=os.path.join(cwd, "checkpoint.swml"), checkpoint=os.path.join(cwd, "checkpoint.npy"))
   ux_initial_from_checkpoint = sw.solution_old.split()[0]
   h_initial_from_checkpoint = sw.solution_old.split()[1]
   
   return ux_initial.vector().array(), h_initial.vector().array(), ux_initial_from_checkpoint.vector().array(), h_initial_from_checkpoint.vector().array()
   
def test_checkpoint():
   
   ux_initial, h_initial, ux_initial_from_checkpoint, h_initial_from_checkpoint = checkpoint()

   print h_initial
   print ux_initial
   print h_initial_from_checkpoint
   print ux_initial_from_checkpoint

   assert(all(h_initial == 0.0))
   assert(all(ux_initial == 0.0))
   assert(all(h_initial_from_checkpoint - 1.0 < 1.0e-14)) # h should be initialised to unity everywhere (from the checkpoint file).
   assert(all(ux_initial_from_checkpoint - 0.0 < 1.0e-7)) # ux should be initialised to (almost) zero everywhere (from the checkpoint file).

if __name__ == '__main__':
   pytest.main(os.path.abspath(__file__))
