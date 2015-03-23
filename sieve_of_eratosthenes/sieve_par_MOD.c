/*
 * Project: COMP 401 - Project03
 *
 * Author: Devin Delfino
 *
 * File Name: sieve_par_MOD.c
 *
 * File Contents: Implementation of a parallel algorithm for the Sieve of Eratosthenes with various modifications for
 *                better efficiency. This implementation takes half of the storage of the original parallel algorithm by
 * 				  removing all even numbers from the sieve. Also, it removes the broadcasting step and reorders the sieve
 *                loops for speed.
 *				         
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
// #define DONT_PRINT

int main(int argc, char* argv[]) {
	int procs;		// the number of processes
	int rank;		// the rank of a given process
	long range = 100;		// the upper bound of the range to find primes (2 to range)
	long sieve_size;
	double elapsed_time;
	int local_prime_count;
	long global_prime_count = 0; // start with 1 to account for the prime '2'

	MPI_Init(&argc, &argv);		// Initialize MPI, breaks into child processes -----------------------------------

	if(argc > 1) {
		range = atoi(argv[1]);
	}

	// declare helper functions so they are accessible by each process
	long get_block_lowest(int,int,long);
	long get_block_highest(int,int,long);
	long get_block_size(int,int,long);
	// long get_block_owner(int,int,long);

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);		// determines the rank of the current process of the communicator
	MPI_Comm_size(MPI_COMM_WORLD, &procs);		// determine the total number of processes in the communicator

	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = MPI_Wtime();

	// ================ MOD ================
	// removing all of the even numbers
	sieve_size = ((range-1) / 2) + 1;
	// ================ MOD ================

	// ================ MOD ================
	long first = get_block_lowest(rank, procs, sieve_size) * 2 + 1;	// gets the first number of the block
	long last = get_block_highest(rank, procs, sieve_size) * 2 + 1;	// gets the last number of the block
	long size = get_block_size(rank, procs, sieve_size);		// gets the size of the block
	// ================ MOD ================

	printf("(Process %d) %ld - %ld, size %ld\n", rank, first, last, size);
	fflush(stdout);

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

	// ================ MOD ================
	// sequential algorithm finding primes between 3 and sqrt(n)
	// given index, actual number = 2 * index + 3
	// given actual number, index = actual number - 3 / 2
	long sqrt_n = (int) sqrt( (double) range);
	char* prime_list;

	prime_list = (char *) malloc (sqrt_n - 2);
	if(prime_list == NULL) {
		printf("Cannot allocate enough memory.\n");
		exit(1);
	}

	for(i=0; i < sqrt_n-2; i++) {
		prime_list[i] = 0;
	}

	int current_prime = 3;
	int current_squared = 9;
	int mark, it;
	while(current_squared <= sqrt_n) {
		for(mark = (current_squared-3)/2; mark < sqrt_n-2; mark += current_prime) {
			prime_list[mark] = 1;
		}

		for(it = ((current_prime-3)/2)+1; it < sqrt_n-2; it++) {
			if(prime_list[it] == 0) {
				current_prime = (it*2)+3;
				break;
			}
		}
		current_squared = current_prime * current_prime;
	}

	// ================ MOD ================

	// ================ MOD ================
	// reorganizing inner and outer loops to improve cache hit rate
	// given index, actual number = first + index * 2
	// given actual number, index = 
	long number;
	long start = 0;
	if(rank == 0) {
		start = 1;
	}

	// iterate through each number in the block
	for(mark = start; mark < size; mark++) {
		number = first + (mark * 2);
		current_prime = 3;
		current_squared = 9;
		while(current_squared <= range+1) {

			if((number % current_prime == 0) && (number >= current_squared)) { // if number at index mark is a multiple of current prime and greater than the square of the current prime
				block[mark] = 1; // mark as composite and break out of while loop to move on to next element in block
				break;
			}
			
			// if number is not composite, find next prime
			for(it = ((current_prime-3)/2)+1; it < sqrt_n-2; it++) {
				if(prime_list[it] == 0) {
					current_prime = (it*2)+3;
					break;
				}
			}
			current_squared = current_prime * current_prime;
		}
	}
	// ================ MOD ================

	local_prime_count = 0;

	if(rank == 0) {
		local_prime_count++; // for the prime '2'
		#ifndef DONT_PRINT
			printf("(%d) 2\n", rank);
			fflush(stdout);
		#endif
	}

	for(it = 1; it < size; it++) {
		if(block[it] == 0 ) {
			local_prime_count++;
			#ifndef DONT_PRINT
				printf("(%d) %ld\n", rank, first + (it*2));
				fflush(stdout);
			#endif
		}
	}

	// Reduction function to increase the global count of prime numbers stored in process 0
	MPI_Reduce(&local_prime_count, &global_prime_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	elapsed_time = MPI_Wtime() - elapsed_time;

	MPI_Finalize();		// Terminates the MPI

	if(rank==0) {
		printf("----------------------------------------\n");
		printf("Total primes between 2 and %ld: %ld\n", range, global_prime_count);
		printf("Time Elapsed: %f seconds\n\n", elapsed_time);
	}

	return 0;
}

// ------ Helper Functions ------------------------------------------
long get_block_lowest(int rank, int procs, long range) {
// Gets the integer that is represented by the first (lowest) index of the block with ID 'rank'
// 	 Parameters: rank - an integer representing the rank of the process
//               procs - an integer representing the number of processes
//               range - an integer representing the maximum number of the sieve
//   Returns: integer that is represented by the lowest index of the block
	return (rank * range)/ procs;
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