// Circuit Satisfiability

#include <mpi.h>
#include <stdio.h>

int main (int argc, char *argv[]) {

	double elapsed_time;
	int solutions;
	int globalsolutions;
	int i;
	int id;						// process rank
	int p; 						// number of processes

	// MPI_Initialized(); 			// checks to see if MPI has been initialized
	// MPI_Init(&argc, &argv);		// initializes MPI
	// MPI_Initialized();

	int check_circuit(int);

	MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = MPI_Wtime();

	solutions = 0;
	for(i = 0; i < 65536; i++) {
		solutions += check_circuit(i);
	}
	elapsed_time = MPI_Wtime() - elapsed_time;
	MPI_Finalize();
	printf("Solutions: %d\n", solutions);
	printf("Time Elapsed: %f seconds\n", elapsed_time);
	return 0;
}

#define EXTRACT_BIT(n,i) ((n& (1<<i))?1:0)

int check_circuit(int z) {
	int v[16];
	int i;

	for(i = 0; i < 16; i++) {
		v[i] = EXTRACT_BIT(z,i);
	}

	if((v[0] || v[1]) && (!v[1] || !v[3]) && (v[2] || v[3])
		&& (!v[3] || !v[4]) && (v[4] || !v[5])
		&& (v[5] || !v[6]) && (v[5] || v[6])
		&& (v[6] || !v[15]) && (v[7] || !v[8])
		&& (!v[7] || !v[13]) && (v[8] || v[9])
		&& (v[8] || !v[9]) && (!v[9] || !v[10])
		&& (v[9] || v[11]) && (v[10] || v[11])
		&& (v[12] || v[13]) && (v[13] || !v[14])
		&& (v[14] || v[15])) {
		// printf("%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d\n",
			// v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
			// v[10],v[11],v[12],v[13],v[14], v[15]);
		// fflush(stdout);
		return 1;
	}
	else {
		return 0;
	}
}