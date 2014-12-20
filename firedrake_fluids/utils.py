import mpi4py

def rank0():
   """ Determine whether or not the calling process has MPI rank 0. """
   if(mpi4py.MPI.COMM_WORLD.Get_rank() == 0):
      return True
   else:
      return False
