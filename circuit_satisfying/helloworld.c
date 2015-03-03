#include <stdio.h>
#include <string.h>
#include <mpi.h>
 
int
main(int argc, char *argv[])
{
    const char *s = "HELLO FROM THE MASTER PROCESS!";
    int par_rank, par_size;
    FILE *fp = NULL;
    char msgbuf[100], filename[100];
 
    /* Init MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &par_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &par_size);
 
    msgbuf[0] = '\0';
 
    /* Broadcast message from master to all other processors. */
    if(par_rank == 0)
    {
        MPI_Bcast((void*)s, strlen(s)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
        strcpy(msgbuf, s);
    }
    else
        MPI_Bcast((void*)msgbuf, strlen(s)+1, MPI_CHAR, 0, MPI_COMM_WORLD);
 
    /* Write the message from the master to a file. */
    sprintf(filename, "%s.%04d.log", argv[0], par_rank);
    if((fp = fopen(filename, "wt")) != NULL)
    {
        fprintf(fp, "Running %s with %d processors.\n", argv[0], par_size);
        fprintf(fp, "This is the log for processor %d.\n", par_rank);
        fprintf(fp, "Message: \"%s\"\n", msgbuf);
        fclose(fp);
    }
 
    /* Finalize MPI */
    MPI_Finalize();
 
    return 0;
}