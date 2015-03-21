/*
 * Project: COMP 401 - Project03
 *
 * Author: Devin Delfino
 *
 * File Name: sieve_seq.c
 *
 * File Contents: Implementation of the sequential algorithm for the Sieve of Eratosthenes.
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>


int main(int argc, char* argv[]) {
	int upper_bound = 100;
	int prime_count = 0;
	clock_t start, stop;
 
	if(argc > 1) {
		upper_bound = atoi(argv[1]);
		// printf("%d\n",atoi(argv[1]));
	}
	
	start = clock();

	char* prime_list;
	prime_list = (char *) malloc (upper_bound-1);
	if(prime_list == NULL) {
		printf("Cannot allocate enough memory.\n");
		exit(1);
	}
	int i;
	for(i=0; i < upper_bound-1; i++) {
		prime_list[i] = 0;
	}

	int current_prime = 2;
	int current_squared = 4;
	int mark, it;
	while(current_squared <= upper_bound) {
		for(mark = current_squared-2; mark < upper_bound-1; mark += current_prime) {
			prime_list[mark] = 1;
		}

		for(it = current_prime-2+1; current_prime < upper_bound-1; it++) {
			if(prime_list[it] == 0) {
				current_prime = it+2;
				break;
			}
		}
		current_squared = current_prime * current_prime;
	}

	printf("Prime Numbers up to %d:\n", upper_bound);
	for(it=0; it < upper_bound-1; it++){
		if(prime_list[it] == 0) {
			prime_count++;
			// printf("%d\n", it+2);
		}	
	}

	stop = clock();


	printf("Total primes between 2 and %d: %d\n", upper_bound, prime_count);
	printf("Time Elapsed: %f seconds\n\n", (double) (stop - start) / CLOCKS_PER_SEC);

	return 1;
}