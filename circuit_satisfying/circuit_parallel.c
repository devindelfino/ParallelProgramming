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
	MPI_Init(&argc, &argv);		// initializes MPI
	// MPI_Initialized();

	int check_circuit(int, int);

	MPI_Comm_rank(MPI_COMM_WORLD, &id);	// determine the rank within a communicator
	MPI_Comm_size(MPI_COMM_WORLD, &p);	// determine the total number of processes in a communicator

	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = MPI_Wtime();

	solutions = 0;
	for(i = id; i < 65536; i += p) {
		solutions += check_circuit(id, i);
	}

	MPI_Reduce(&solutions, &globalsolutions, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	elapsed_time = MPI_Wtime() - elapsed_time;

	printf("Process %d is done\n", id);
	fflush(stdout);
	MPI_Finalize(); 		// shuts off MPI

	if(id==0) {
		printf("Solutions: %d\n", globalsolutions);
	}
	printf("Time Elapsed: %f seconds", elapsed_time);
	return 0;
}

#define EXTRACT_BIT(n,i) ((n& (1<<i))?1:0)

int check_circuit(int id, int z) {
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
		printf("%d) %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d\n", id,
			v[0],v[1],v[2],v[3],v[4],v[5],v[6],v[7],v[8],v[9],
			v[10],v[11],v[12],v[13],v[14], v[15]);
		fflush(stdout);
		return 1;
	}
	else {
		return 0;
	}
}