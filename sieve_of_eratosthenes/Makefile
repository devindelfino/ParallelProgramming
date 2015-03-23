#targets-----------------------------
sequential_sieve:	
				gcc sieve_seq.c -o sequential_sieve -lm

parallel_sieve:	
				mpicc sieve_par.c -o parallel_sieve -lm

mod_sieve:	
				mpicc sieve_par_MOD.c -o mod_sieve -lm
				
clean:
	rm *_sieve