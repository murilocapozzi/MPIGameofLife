/*

Jogo da Vida em MPI

Autor: Murilo Capozzi dos Santos
RA: 149425
Turma: I

*/

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define BOARD_SIZE 2048
#define GENERATIONS 2000

void initGrid(float **grid, int rank, int size){
    int i, j;

    for(i = rank; i < BOARD_SIZE; i += size){
        grid[i] = malloc(BOARD_SIZE*sizeof(float));
        for(j = 0; j < BOARD_SIZE; j++){
            if(i == 1 && j == 2 ||
               i == 2 && j == 3 ||
               i == 3 && (j == 1 || j == 2 || j == 3) ||
               i == 10 && (j == 31 || j == 32) ||
               i == 11 && (j == 30 || j == 31) ||
               i == 12 && j == 31)
               grid[i][j] = 1.0;
            else grid[i][j] = 0.0;
        }
    }
}

int getCoord(int n){

    while(n % BOARD_SIZE < 0)
        n += BOARD_SIZE;

    return n % BOARD_SIZE;
}

int getNeighbors(float *grid, float *recv_up, float *recv_down, int j){
    int k, quantity = 0;

    for(k = -1; k <= 1; k++){
        if(recv_up[getCoord(j + k)] > 0.0)
            quantity++;
        if(recv_down[getCoord(j + k)] > 0.0)
            quantity++;
        if(k != 0 && grid[getCoord(j + k)] > 0.0)
            quantity++;
    }

    return quantity;
}


float meanNeighbors(float *grid, float *recv_up, float *recv_down, int j){
    int k;
    float sum = 0.0;

    for(k = -1; k <= 1; k++){
        
        if(k != 0 && grid[getCoord(j + k)] > 0.0)
            sum += grid[getCoord(j + k)];
        sum += recv_up[getCoord(j + k)] + recv_down[getCoord(j + k)];
    }
    
    return (sum / 8.0);
}


void gameOfLife(float **grid, float **new_grid, int rank, int size){
    int i, j, quantityNeighbors, quantityAlive = 0, totalAlive, up_want_i, down_want_i;
    int neighbor_up = (rank - 1 + size) % size;
    int neighbor_down = (rank + 1) % size;
    float value, recv_up[BOARD_SIZE], recv_down[BOARD_SIZE];

    for(i = rank; i < BOARD_SIZE; i += size){

        MPI_Sendrecv(&i, 1, MPI_INT, neighbor_down, 0, &up_want_i, 1, MPI_INT, neighbor_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(&i, 1, MPI_INT, neighbor_up, 1, &down_want_i, 1, MPI_INT, neighbor_down, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Barrier(MPI_COMM_WORLD);
        
        MPI_Sendrecv(grid[getCoord(down_want_i - 1)], BOARD_SIZE, MPI_FLOAT, neighbor_down, 0, recv_up, BOARD_SIZE, MPI_FLOAT, neighbor_up, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        MPI_Sendrecv(grid[getCoord(up_want_i + 1)], BOARD_SIZE, MPI_FLOAT, neighbor_up, 1, recv_down, BOARD_SIZE, MPI_FLOAT, neighbor_down, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for(j = 0; j < BOARD_SIZE; j++){
            
            quantityNeighbors = getNeighbors(grid[i], recv_up, recv_down, j);

            value = grid[i][j];

            if(value > 0 && (quantityNeighbors == 2 || quantityNeighbors == 3)){
                new_grid[i][j] = 1.0;
                quantityAlive++;
            }
            else if(value == 0 && quantityNeighbors == 3){
                new_grid[i][j] = meanNeighbors(grid[i], recv_up, recv_down, j);
                quantityAlive++;
            }
            else
                new_grid[i][j] = 0.0;
        }

    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&quantityAlive, &totalAlive, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if(rank == 0) printf("%d vivas\n", totalAlive);
}

void freeGrid(float **grid, int rank, int size){
    int i;
    for(i = rank; i < BOARD_SIZE; i+=size){
        free(grid[i]);
    }
}

int main(int argc, char **argv){

    MPI_Init(&argc, &argv);

    int i, rank, size;
    float **actual, **new, **temp;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    actual = malloc(BOARD_SIZE*sizeof(float*));
    new = malloc(BOARD_SIZE*sizeof(float*));

    initGrid(actual, rank, size);
    initGrid(new, rank, size);

    for(i = 0; i < GENERATIONS; i++){

        if(rank == 0) printf("Geração %d: ", i+1);

        gameOfLife(actual, new, rank, size);
        
        temp = actual;
        actual = new;
        new = temp;
    }

    freeGrid(actual, rank, size);
    freeGrid(new, rank, size);
    free(actual);
    free(new);
    
    MPI_Finalize();

    return 0;
}