all:
	clear

compila:
	mpicc -o mpi_rainbow mpi_rainbow.c -O3

um:
	time mpirun -np 1 ./mpi_rainbow

dois:
	time mpirun -np 2 ./mpi_rainbow

quatro:
	time mpirun -np 4 ./mpi_rainbow

oito:
	time mpirun -np 8 ./mpi_rainbow