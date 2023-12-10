all:
	clear

compila:
	mpicc -o mpi_rainbow mpi_rainbow.c -O3

um:
	mpirun -np 1 ./mpi_rainbow

dois:
	mpirun -np 2 ./mpi_rainbow

quatro:
	mpirun -np 4 ./mpi_rainbow

oito:
	mpirun -np 8 ./mpi_rainbow