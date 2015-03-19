// sieve_par.c

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// function declarations
int get_block_lowest(int,int,int);
int get_block_highest(int,int,int);
int get_block_size(int,int,int);
int get_block_owner(int,int,int);

int main(int argc, char* argv[]) {
	int procs;		// the number of processes
	int rank;		// the rank of a given process
	int range = 100;		// the upper bound of the range to find primes (2 to range)

	MPI_Init(&argc, &argv);		// Initialize MPI, breaks into child processes -----------------------------------

	if(argc > 1) {
		range = atoi(argv[1]);
	}

	// declare helper functions so they are accessible by each process
	int get_block_lowest(int,int,int);
	int get_block_highest(int,int,int);
	int get_block_size(int,int,int);
	int get_block_owner(int,int,int);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);		// determines the rank of the current process of the communicator
	MPI_Comm_size(MPI_COMM_WORLD, &procs);		// determine the total number of processes in the communicator

	int first = get_block_lowest(rank, procs, range);
	int last = get_block_highest(rank, procs, range);
	int size = get_block_size(rank, procs, range);

	MPI_Finalize();		// Terminates the MPI
	return 1;
}


// ------ Helper Functions ------------------------------------------
int get_block_lowest(int rank, int procs, int range) {
	return (rank * range)/ procs;
}

int get_block_highest(int rank, int num_procs, int range) {
	return (((rank + 1) * range) / procs) - 1;
}

int get_block_size(int rank, int procs, int range) {
	return get_block_highest(rank,procs,range) - get_block_lowest(rank,procs,range) + 1;
}

int get_block_owner(int ind, int num_procs, int range) {
	return (procs * (ind + 1) - 1) / range;
}