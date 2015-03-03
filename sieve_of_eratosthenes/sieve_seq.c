// sieve_seq.c

#include <stdio.h>
#include <string.h>
#include <stdlib.h>


int main(int argc, char* argv[]) {
	int upper_bound = 100;

	if(argc > 1) {
		upper_bound = atoi(argv[1]);
		// printf("%d\n",atoi(argv[1]));
	}
	
	int prime_list[upper_bound+1];
	memset(prime_list, 0, sizeof(prime_list));

	int current = 2;
	int current_squared = 4;
	while(current_squared <= upper_bound) {
		for(int mark = current_squared; mark < upper_bound+1; mark += current) {
			prime_list[mark] = 1;
		}

		for(int it = current+1; current < upper_bound+1; it++) {
			if(prime_list[it] == 0) {
				current = it;
				break;
			}
		}
		current_squared = current * current;
	}

	printf("Prime Numbers up to %d:\n", upper_bound);
	for(int prime=2; prime < upper_bound+1; prime++){
		if(prime_list[prime] == 0) {
			printf("%d\n", prime);
		}	
	}
	return 1;
}