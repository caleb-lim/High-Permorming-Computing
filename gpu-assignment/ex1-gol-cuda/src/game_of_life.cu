#include <stdio.h>
#include <stdlib.h>
#include "common.h"

#define INCLUDE_CPU_VERSION 
#include "game_of_life.c"

#define USEPNG

#define NTHREADS 1024

#define CUDA_CHECK_ERROR(X)({\
    if((X) != cudaSuccess){\
        fprintf(stderr, "CUDA error (%s:%d): %s\n", __FILE__, __LINE__, cudaGetErrorString((X)));\
        exit(1);\
    }\
})

#define MALLOC_CHECK_ERROR(X)({\
    if ((X) == 0){\
        fprintf(stderr, "Malloc error (%s:%d): %i\n", __FILE__, __LINE__, (X));\
        exit(1);\
    }\
})

__global__ void game_of_life_step(int *current_grid, int *next_grid, int n, int m){
    unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    
    int neighbours;
    int n_i[8], n_j[8];

    if (idx < n * m){
        unsigned int i = idx / m;
        unsigned int j = idx % m;
        
        // count the number of neighbours, clockwise around the current cell.
        neighbours = 0;
        
        n_i[0] = i - 1; n_j[0] = j - 1;
        n_i[1] = i - 1; n_j[1] = j;
        n_i[2] = i - 1; n_j[2] = j + 1;
        n_i[3] = i;     n_j[3] = j + 1;
        n_i[4] = i + 1; n_j[4] = j + 1;
        n_i[5] = i + 1; n_j[5] = j;
        n_i[6] = i + 1; n_j[6] = j - 1;
        n_i[7] = i;     n_j[7] = j - 1;

        if(n_i[0] >= 0 && n_j[0] >= 0 && current_grid[n_i[0] * m + n_j[0]] == ALIVE) neighbours++;
        if(n_i[1] >= 0 && current_grid[n_i[1] * m + n_j[1]] == ALIVE) neighbours++;
        if(n_i[2] >= 0 && n_j[2] < m && current_grid[n_i[2] * m + n_j[2]] == ALIVE) neighbours++;
        if(n_j[3] < m && current_grid[n_i[3] * m + n_j[3]] == ALIVE) neighbours++;
        if(n_i[4] < n && n_j[4] < m && current_grid[n_i[4] * m + n_j[4]] == ALIVE) neighbours++;
        if(n_i[5] < n && current_grid[n_i[5] * m + n_j[5]] == ALIVE) neighbours++;
        if(n_i[6] < n && n_j[6] >= 0 && current_grid[n_i[6] * m + n_j[6]] == ALIVE) neighbours++;
        if(n_j[7] >= 0 && current_grid[n_i[7] * m + n_j[7]] == ALIVE) neighbours++;

        if(current_grid[i*m + j] == ALIVE && (neighbours == 2 || neighbours == 3)){
            next_grid[i*m + j] = ALIVE;
        } else if(current_grid[i*m + j] == DEAD && neighbours == 3){
            next_grid[i*m + j] = ALIVE;
        }else{
            next_grid[i*m + j] = DEAD;
        }

    }

}

int* game_of_life(const int *initial_state, int n, int m, int nsteps){
    // Allocate memory for the grids
    int *grid = (int *) malloc(sizeof(int) * n * m);
    int *updated_grid = (int *) malloc(sizeof(int) * n * m);
    if(!grid || !updated_grid){
        printf("Error while allocating memory.\n");
        exit(1);
    }

    // Copy the initial state to the grid
    memcpy(grid, initial_state, sizeof(int) * n * m);

    // Setup kernel configuration
    // const int BLOCK_SIZE = 256;
    unsigned int nBlocks = (n * m + NTHREADS - 1) / NTHREADS;

    // Allocate device memory for the grids
    int *dev_grid, *dev_updated_grid;
    CUDA_CHECK_ERROR(cudaMalloc(&dev_grid, sizeof(int) * n * m));
    CUDA_CHECK_ERROR(cudaMalloc(&dev_updated_grid, sizeof(int) * n * m));

    // Copy the grid to the device
    CUDA_CHECK_ERROR(cudaMemcpy(dev_grid, grid, sizeof(int) * n * m, cudaMemcpyHostToDevice));

    // Run the simulation for nsteps steps
    for(int step = 0; step < nsteps; step++) {
        // Launch the kernel
        game_of_life_step<<<nBlocks, NTHREADS>>>(dev_grid, dev_updated_grid, n, m);
        CUDA_CHECK_ERROR(cudaGetLastError());

        // Swap the grids
        int *temp = dev_grid;
        dev_grid = dev_updated_grid;
        dev_updated_grid = temp;
    }

    // Copy the final grid back to the host
    CUDA_CHECK_ERROR(cudaMemcpy(grid, dev_grid, sizeof(int) * n * m, cudaMemcpyDeviceToHost));

    // Free device memory
    CUDA_CHECK_ERROR(cudaFree(dev_grid));
    CUDA_CHECK_ERROR(cudaFree(dev_updated_grid));
    free(updated_grid);

    return grid;
}

int main(int argc, char **argv)
{
    printf("GPU version");
    struct Options *opt = (struct Options *) malloc(sizeof(struct Options));
    getinput(argc, argv, opt);
    int n = opt->n, m = opt->m, nsteps = opt->nsteps;
    int *initial_state = (int *) malloc(sizeof(int) * n * m);
    if(!initial_state){
        printf("Error while allocating memory.\n");
        return -1;
    }
    generate_IC(opt->iictype, initial_state, n, m);
    struct timeval start_cuda, start;
    start_cuda = init_time();
    int *final_state_cuda = game_of_life(initial_state, n, m, nsteps);
    float elapsed_cuda = get_elapsed_time(start_cuda);
    printf("Finnished CUDA in %f ms\n", elapsed_cuda);

    start = init_time();
    int *final_state = cpu_game_of_life(initial_state, n, m, nsteps);
    float elapsed = get_elapsed_time(start);
    printf("Finnished GOL in %f ms\n", elapsed);

    compare(final_state, final_state_cuda, n, m);
    visualise(VISUAL_ASCII, 100, final_state, n, m);
    visualise(VISUAL_ASCII, 100, final_state_cuda, n, m);

    free(final_state_cuda);
    free(final_state);
    free(initial_state);
    free(opt);
    return 0;
}
