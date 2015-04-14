/*
 * Project: COMP 401 - Project03
 *
 * Author: Devin Delfino
 *
 * File Name: sieve_par.c
 *
 * File Contents: Implementation of a parallel algorithm for the Sieve of Eratosthenes using MPI for the
 *                message passing system. This is a direct parallelized version of the sequential algorithm.         
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define DONT_PRINT

int main(int argc, char* argv[]) {
	int procs;						// the number of processes
	int rank;						// the rank of a given process
	long range = 100;				// the upper bound of the range to find primes (2 to range)
	double elapsed_time;			// timekeeping
	int local_prime_count;			// local prime count in a particular process
	long global_prime_count = 0;	// total prime count 

	MPI_Init(&argc, &argv);		// Initialize MPI, breaks into child processes -----------------------------------

	if(argc > 1) {
		range = atol(argv[1]); }

	// ------------------------------------------------------------------
	// declare helper functions so they are accessible by each process
	long get_block_lowest(int,int,long);
	long get_block_highest(int,int,long);
	long get_block_size(int,int,long);
	// long get_block_owner(int,int,long);
	// ------------------------------------------------------------------

	// ------------------------------------------------------------------
	// MPI_Comm_rank()
	// determines the rank of the current process
	// Parameters:
	// MPI_COMM_WORLD - object that provides the environment for message passing among processes
	//          &rank - pass by reference, updates 'rank' to the rank of the current process
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	// MPI_Comm_size()
	// determines the total number of processes
	// Parameters:
	// MPI_COMM_WORLD - object that provides the environment for message passing among processes
	//          &procs - pass by reference, updates 'size' to the total number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	// MPI_Barrier()
	// no process can proceed beyond a barrier until all processes have reached it
	// ensures all processes enter the same block of code at the same time 
	// Parameters:
	// MPI_COMM_WORLD - object that provides the environment for message passing among processes
	MPI_Barrier(MPI_COMM_WORLD);

	// MPI_Wtime()
	// returns number of elapsed seconds
	elapsed_time = MPI_Wtime();
	// ------------------------------------------------------------------

	// ------------------------------------------------------------------
	// checks for too many processes - this will result in a higher need for communication, decreasing efficiency
	// ensures n/p > sqrt(n)
	if(rank == 0) {
		int proc0_size = (range - 1)/procs;	// size of the 0th process
		if( (proc0_size) < (int) sqrt( (double) range)) {
			printf("There are currently %d processes.\n", procs);
			printf("If n = %ld, there must be less than %d processes in order to minimize communication.\n", range, (int)sqrt((double)range));
		}
		MPI_Finalize(); // exit MPI
		exit(0);
	}
	// ------------------------------------------------------------------

	// ------------------------------------------------------------------
	// parallelizing step 1 of sequential algorithm 
	// 1. Create a list of natural numbers 2, 3, 4, ... , n, all marked 0

	long first = get_block_lowest(rank, procs, range - 1) + 2;	// gets the first number represented by the block
	long last = get_block_highest(rank, procs, range - 1) + 2;	// gets the last number represented by the block
	long size = get_block_size(rank, procs, range - 1);			// gets the size of the block

	printf("(Process %d) %ld - %ld, size %ld\n", rank, first, last, size);
	fflush(stdout);

	// Creates a block of the entire array to be handled by the current process
	char* block;
	block = (char *) malloc (size);
	if(block == NULL) {
		if(rank == 0) {
			printf("Cannot allocate enough memory.\n");
		}
		MPI_Finalize();
		exit(0);
	}
	long i;
	for(i=0; i < size; i++) {
		block[i] = 0;
	}
	// ------------------------------------------------------------------

	// ------------------------------------------------------------------
	// parallelizing step 2 of sequential algorithm
	// 2. Set k = 2, the first unmarked number on the list 
	int current_prime = 2;
	int current_squared = 4;
	// ------------------------------------------------------------------
	
	// ------------------------------------------------------------------
	// parallelizing step 3a of sequential algorithm ----------------------------------
	// 3a. Mark all multiples of k between k^2 and n

	// Repeat until k^2 > n
	int multiple;	// the local index indicating the first multiple of the current prime in the block
	int rem;		// the remainder of first % current_prime
	int mark, it;   // iterators
	while(current_squared <= range+1) {
		// finding the first multiple of k in this block (greater than or equal to k^2)
		if(current_squared <= first) { // k^2 already occurred, so we just mark the first multiple of the current prime
			rem = first % current_prime;
			if(rem == 0) {
				multiple = 0;
			}
			else { // rem > 0
				multiple = current_prime - rem;
			}
		}
		else { // current_squared > first, meaning we start marking at k^2
			multiple = current_squared - first;
		}

		// set the first mark to the first multiple of the current prime and increase by the current prime
		
		for(mark = multiple; mark < size; mark += current_prime) {
			block[mark] = 1;
		}


		// ------------------------------------------------------------------
		// parallelizing step 3b of sequential algorithm ----------------------------------
		// 3b. Set k to the smallest unmarked number greater than the current k
		if(rank == 0) {
			for(it = current_prime-2+1; it < size; it++) {
				if(block[it] == 0) {
					current_prime = it + 2;
					break;
				}
			}
		}

		// MPI_Bcast()
		// allows a process to broadast a data item to all other processes in the communicator (local communication)
		// Parameters:
		//  	    void *buffer - address of first element to be broadcasted
		//             int count - number of elements to be broadcasted
		// 	   MPI_Datatype type - type of element
		// 			    int root - rank of process performing broadcast
		// MPI_Comm communicator - communicator object
		MPI_Bcast(&current_prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
		
		current_squared = current_prime * current_prime;
		// ------------------------------------------------------------------

	}
	// ------------------------------------------------------------------

	// ------------------------------------------------------------------
	// Iterates through the block, counts the unmarked values, 
	// and increases the local prime count for the process
	local_prime_count = 0;
	for(it = 0; it < size; it++) {
		if(block[it] == 0 ) {
			local_prime_count++;
			#ifndef DONT_PRINT
				printf("(%d) %ld\n", rank, first + it);
				fflush(stdout);
			#endif
		}
	}
	// ------------------------------------------------------------------

	// ------------------------------------------------------------------
	// Reduction function to increase the global count of 
	// prime numbers stored in process 0

	// MPI_Reduce()
	// performs reduction operation on values submitted by all processes in a communicator (global communication)
	// Parameters:
	//  	   void *operand - address of first reduction element
	// 			void *result - address of first reduction result
	//             int count - number of reductions to perform
	// 	   MPI_Datatype type - type of element
	// 	     MPI_Op operator - reduction operator
	// 			    int root - process receiving result
	// MPI_Comm communicator - communicator object
	MPI_Reduce(&local_prime_count, &global_prime_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	elapsed_time = MPI_Wtime() - elapsed_time;

	MPI_Finalize();		// Terminates the MPI
	// ------------------------------------------------------------------

	// ------------------------------------------------------------------
	// Print number of primes and elapsed time
	if(rank==0) {
		printf("----------------------------------------\n");
		printf("Total primes between 2 and %ld: %ld\n", range, global_prime_count);
		printf("Time Elapsed: %f seconds\n\n", elapsed_time);
	}
	// ------------------------------------------------------------------

	return 0;
}


// ------ Helper Functions ------------------------------------------
long get_block_lowest(int rank, int procs, long range) {
// Gets the integer that is represented by the first (lowest) index of the block with ID 'rank'
// 	 Parameters: rank - an integer representing the rank of the process
//               procs - an integer representing the number of processes
//               range - an integer representing the maximum number of the sieve
//   Returns: integer that is represented by the lowest index of the block
	return (rank * range) / procs;
}

long get_block_highest(int rank, int procs, long range) {
// Gets the integer that is represented by the last (highest) index of the block with ID 'rank'
// 	 Parameters: rank - an integer representing the rank of the process
//               procs - an integer representing the number of processes
//               range - an integer representing the maximum number of the sieve
//   Returns: integer that is represented by the highest index of the block
	return (((rank + 1) * range) / procs) - 1;
}

long get_block_size(int rank, int procs, long range) {
// Gets the size of the block with ID 'rank'
// 	 Parameters: rank - an integer representing the rank of the process
//               procs - an integer representing the number of processes
//               range - an integer representing the maximum number of the sieve
//   Returns: integer that is represented by the lowest index of the block
	return get_block_highest(rank,procs,range) - get_block_lowest(rank,procs,range) + 1;
}

// int get_block_owner(int ind, int procs, long range) {
// 	return (procs * (ind + 1) - 1) / range;
// }
