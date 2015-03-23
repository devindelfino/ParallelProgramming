Project 03: The Sieve of Eratosthenes
=====================================
*Devin Delfino*

Description
-----------
[The Sieve of Eratosthenes](http://en.wikipedia.org/wiki/Sieve_of_Eratosthenes) is an algorithm that computes all of the prime number up to and including a given limit. This is accomplished by marking every multiple of each prime (starting with 2) as composite. For each prime *p*, the algorithm marks each multiple of *p* starting with *p^2*, and then does the same for each succesive *p* until *p^2* is greater than the given limit.

The parallel algorithm uses a Block Decomposition scheme to break the large of number between 2 and the maximum limit into contiguous blocks, where each block is handled by a particular process. This parallel program uses the [Message Passing Interface (MPI)](http://en.wikipedia.org/wiki/Message_Passing_Interface) system to handle the parallel aspects of the program.

Documentation
-------------

#####Installing Dependencies
In order to compile code using MPI, you must do the following:

1. [Download](http://www.open-mpi.org/software/ompi/v1.8/) Open MPI
2. [Install](https://wiki.helsinki.fi/display/HUGG/Installing+Open+MPI+on+Mac+OS+X) Open MPI (the link is for Mac OS X)


#####Running the Sequential Program

To compile code...

	make sequential_sieve

To run executable with no parameters...
		     
	./sequential_sieve

*Defaults to:*
*Maximum number of the sieve = 100 (finds primes between 2 and 100, inclusive)*

To run executable with parameters...
		     
	./sequential_sieve N

*Maximum number of the sieve = __N__ (finds primes between 2 and n, inclusive)*

To remove all executable files (required before recompiling code after changes are made)

	make clean

#####Running the Parallel Program

To compile code...

	make parallel_sieve

To run executable with no parameters...

	mpirun parallel_sieve 

*Defaults to:*
*Number of processes = Number of CPU cores*
*Maximum number of the sieve = 100 (finds primes between 2 and 100)*

To run executable with parameters...
		     
	mpirun -np P parallel_sieve N

*Number of processes = __P__*
*Maximum number of the sieve = __N__ (finds primes between 2 and n, inclusive)*

To remove all of the executable files (required before recompiling code after changes are made)

	make clean

Complexity Analysis
-------------------
__Sequential Algorithm__
O(n log( log( n ) ) )

__Parallel Algorithm__
O( X (n log log n) / p + (sqrt(n)/log sqrt(n)) Y ceil(log p) )

Contents
--------
The *sieve_of_eratosthenes* directory of the *ParallelProgramming* repository contains the following files:

* *sieve_seq.c* - A sequential implementation of the sieve of eratosthenes algorithm
* *sieve_par.c* - A parallel implementation of the sieve of eratosthenes algorithm

References
----------
1. Parallel Programming in C with MPI and OpenMP by Michael J. Quinn


Notes
-----

Parallel Algorithm Design (Reference 1, Chapter 3)
--------------------------------------------------
* The task/channel model represetns a parallel computation as a set of tass that may interact with each other by sending messages through channels
* A task is a program, its local memory, and a collection of I/O ports
* A channel is a message queue that connects one task's output port with another task's input port
* Binomial trees are the most common communication pattern in parallel algorithm design

Foster's Design Methodology
-------------------

1. Partitioning

* def: the porcess of dividing the computation and the data into pieces
* domain decomposition is the parallel algorithm design approach in which we first divide the data into pieces and then determine how to associate computations with the data (it is best to maximize primitive tasks)
* functional decomposition is the strategy in which we first divide the computation into pieces nd then determine how to associate data items with the individual computations
* checklist: Page 67
 
2. Communication

* local communication - when a task needs values from a small number of other tasks in order to perform a computation, we create channels from the tasks supplying data to the task consuming the data
* global communication - when a significant number of the primitive tasks must contribute data in order to perform a computation
* Minimize parallel overhead is very important
* Checklist: commun. operations are balanced among tasks, each task communicates iwth only a small number of neighbors, task can perform their communications and computation concurrently

3. Agglomeration

* def: the process of gropuing tasks into larger tasks in order to improve perforamce or simplify programming. 
* agglomeration lowers communication overhead
* increasing locality of the parallel algorithm is agglomerateing primitive tasks that communicate with each other, resulting in communication between them is eliminated
* can also gropu together sending and receiving tasks, reducing number of messages being sent
* secondary goal is to maintain scalability of the design so it can be later applied to a computer with more processors
* checklist: pg 69

4. Mapping

* def: process of assigning tasks to processors
* Processor utilization is the average percentage of time the system's procesors are actively executing tasks necessary for the solution of the problem

Message-Passing Programming (Reference 1, Chapter 4)
--------------------------------------------------
* rank of a process is its position in the overall order
MPI Functions
-------------
* MPI_Init() - Initialize MPI
* MPI_Comm_rank() - determines a process's ID number
* MPI_Comm_size() - finds the number of processes
* MPI_Reduce() - performs a reduction operation
* MPI_Finalize() - shuts down MPI
* MPI_Barrier() - performs a barrier synchronization
* MPI_Wtime() - determines the time
* MPI_Wtick() - finds the accuracy of the timer

The Sieve of Eratosthenes (Reference 1, Chapter 5)
--------------------------------------------------
Data Decomposition Options

####Interleaved

* process 0 handles 2, 2 + p, 2 + 2p, ..., process 1 handles 3, 3 + p, 3 + 2p, ...
* Disadvantage - Unbalanced work between processes, requires some reduction operation

####Block

* Key Question 1: What is the range of elements controlled y a particular process?
* Key Question 2: Which process controls a particular element?
* If the array has n elements and there are p available processes, assign each process a block of floor(n/p) or ceiling(n/p) elements ( if n is divisible by p, then each process gets assigned n/p elements)

######Block Allocation Method

1. The first element controlled by process i is floor( (i*n)/p ) 
2. The final element controlled by process i is floor( ((i+1)*n)/p ) - 1
3. The process controlling a particular array element j is floor( (p(j+1)-1)/n )