build:
	mpicc -fopenmp main_multi_complete.c -o main_multi_omp
clean:
	rm main_multi_omp

