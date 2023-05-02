def get_MPI_data(comm=None, size=None, rank=None):
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
        size = comm.Get_size() # total number of processors
        rank = comm.Get_rank() # the no. of current processor
        # print('multi-threading, processor %i out of %i' % (rank, size))
        if size > 1:
            MPI_enabled = True
        elif size == 1:
            MPI_enabled = False
        else: # size <= 0
            raise ValueError("Total number of processors should not be less than 1.")
    except ImportError:
        MPI_enabled = False
    return MPI_enabled, comm, size, rank