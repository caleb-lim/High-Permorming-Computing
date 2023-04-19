#include <stdio.h>
#include <stdlib.h>
#include "common.h"

#define INCLUDE_CPU_VERSION 
#include "game_of_life.c"

#define USEPNG

#define NTHREADS 1024
#define BLOCK_SIZE 256
#define HALO_SIZE 1

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


// Split the grid into subgrids with each subgrid having a halo
int* split_grid_to_subgrids(int** grid, int n, int m){
    int num_blocks_x = (m + 2 * HALO_SIZE + BLOCK_SIZE - 1) / BLOCK_SIZE;
    int num_blocks_y = (n + 2 * HALO_SIZE + BLOCK_SIZE - 1) / BLOCK_SIZE;
    
    int block_index = 0;
    int num_subgrid_blocks = num_blocks_x * num_blocks_y;
    int ***blocks = malloc(num_subgrid_blocks * sizeof(int **));


    for (int block_y = 0; block_y < num_blocks_y; block_y++) {
        for (int block_x = 0; block_x < num_blocks_x; block_x++) {
            // Calculate the block's start and end indices
            int block_start_x = block_x * BLOCK_SIZE - HALO_SIZE;
            int block_start_y = block_y * BLOCK_SIZE - HALO_SIZE;
            int block_end_x = (block_x + 1) * BLOCK_SIZE + HALO_SIZE - 1;
            int block_end_y = (block_y + 1) * BLOCK_SIZE + HALO_SIZE - 1;


            // Adjust the indices to handle edge cases
            if (block_start_x < 0) {
                block_start_x = 0;
            }
            if (block_start_y < 0) {
                block_start_y = 0;
            }
            if (block_end_x >= m) {
                block_end_x = m - 1;
            }
            if (block_end_y >= n) {
                block_end_y = n - 1;
            }

            // Calculate the size of the block
            // int block_width = block_end_x - block_start_x + 1;
            // int block_height = block_end_y - block_start_y + 1;

            //print all varibale above for debuugging
            printf("block_start_x: %d, block_start_y: %d, block_end_x: %d, block_end_y: %d\n", block_start_x, block_start_y, block_end_x, block_end_y);

            // Allocate memory for the block
            int **block = malloc(block_height * sizeof(int *));
            for (int i = 0; i < block_height; i++) {
                block[i] = malloc(block_width * sizeof(int));
            }

            // Copy the block from the grid
            for (int i = block_start_y; i <= block_end_y; i++) {
                for (int j = block_start_x; j <= block_end_x; j++) {
                    block[i - block_start_y][j - block_start_x] = grid[i][j];
                }
            }

            // Copy the block's data into the block array
            for (int i = 0; i < block_height; i++) {
                for (int j = 0; j < block_width; j++) {
                        block[i][j] = grid[block_start_y + i][block_start_x + j];      
                  }
            }
           
            // Insert the block into the subgrid blocks array
            blocks[block_y][block_x] = block;
            block_index++;
            

            //Free the block
            
            for (int i = 0; i < block_height; i++) {
                free(block[i]);
            }
            free(block);
            
        }
    return subgrids;
}

__global__ void game_of_life_step(int *current_grid, int *next_grid, int n, int m){
    
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if(i < n && j < m){
        int alive_neighbours = 0;
        for(int k = -1; k <= 1; k++){
            for(int l = -1; l <= 1; l++){
                if(k == 0 && l == 0) continue;
                if(i + k >= 0 && i + k < n && j + l >= 0 && j + l < m){
                    alive_neighbours += current_grid[(i + k) * m + (j + l)];
                }
            }
        }

        if(current_grid[i * m + j] == 1){
            if(alive_neighbours < 2 || alive_neighbours > 3){
                next_grid[i * m + j] = 0;
            }else{
                next_grid[i * m + j] = 1;
            }
        }else{
            if(alive_neighbours == 3){
                next_grid[i * m + j] = 1;
            }else{
                next_grid[i * m + j] = 0;
            }
        }
    }
}

// int* game_of_life(const int **initial_state, int n, int m, int nsteps){
//     int *current_grid, *next_grid;
//     int *d_current_grid, *d_next_grid;
//     int *final_state;

//     int size = n * m * sizeof(int);
//     current_grid = (int *) malloc(size);
//     next_grid = (int *) malloc(size);
//     final_state = (int *) malloc(size);

//     memcpy(current_grid, initial_state, sizeof(int) * n * m);


//     MALLOC_CHECK_ERROR(current_grid);
//     MALLOC_CHECK_ERROR(next_grid);
//     MALLOC_CHECK_ERROR(final_state);
    
//     // Allocate device memory for the grids
//     CUDA_CHECK_ERROR(cudaMalloc((void **) &d_current_grid, size));
//     CUDA_CHECK_ERROR(cudaMalloc((void **) &d_next_grid, size));
    
//     // Copy the grid to the device
//     CUDA_CHECK_ERROR(cudaMemcpy(d_current_grid, current_grid, size, cudaMemcpyHostToDevice));

//     dim3 dimBlock(NTHREADS, NTHREADS);
//     dim3 dimGrid((n + dimBlock.x - 1) / dimBlock.x, (m + dimBlock.y - 1) / dimBlock.y);
//     int blocks_per_dim_x = (n + BLOCK_SIZE - 1)/BLOCK_SIZE;
//     int blocks_per_dim_y  = (m + BLOCK_SIZE - 1)/BLOCK_SIZE;



//     for(int step = 0; step < nsteps; step++){
       
//         //Split the grids into subgrids with each subgrid having a "halo" and compute the next state of each subgrid
//         for(int i = 0; i < blocks_per_dim_x; i++){
//             for(int j = 0; j < blocks_per_dim_y; j++){    
//                 dim3 dimBlock(BLOCK_SIZE, BLOCK_SIZE);
//                 dim3 dimGrid(1, 1);
//                 game_of_life_step<<<dimGrid, dimBlock>>>(d_current_grid, d_next_grid, n, m);
//             }
//         }

//         int *tmp = d_current_grid;
//         d_current_grid = d_next_grid;
//         d_next_grid = tmp;
//     }

//     CUDA_CHECK_ERROR(cudaMemcpy(final_state, d_current_grid, size, cudaMemcpyDeviceToHost));

//     CUDA_CHECK_ERROR(cudaFree(d_current_grid));
//     CUDA_CHECK_ERROR(cudaFree(d_next_grid));

//     free(current_grid);
//     free(next_grid);

//     return final_state;
// }

int main(int argc, char **argv)
{
    printf("GPU version");
    struct Options *opt = (struct Options *) malloc(sizeof(struct Options));
    getinput(argc, argv, opt);
    int n = opt->n, m = opt->m;
    int *initial_state = (int *) malloc(sizeof(int) * n * m);

    int **game_of_life_grid = (int**)malloc(n * sizeof(int*));
    for(int i = 0; i < n; i++){
        game_of_life_grid[i] = (int*)malloc(m * sizeof(int));
    }

    if(!initial_state){
        printf("Error while allocating memory.\n");
        return -1;
    }
    
    generate_IC(opt->iictype, initial_state, n, m);
    convert2D(game_of_life_grid, initial_state, n, m); 
    
    int* subgrids = split_grid_to_subgrids(game_of_life_grid, n, m);
    free(subgrids);
    // struct timeval start_cuda, start;
    // start_cuda = init_time();
    // int **final_state_cuda = game_of_life(game_of_life_grid, n, m, nsteps);

    // float elapsed_cuda = get_elapsed_time(start_cuda);
    // printf("Finnished CUDA in %f ms\n", elapsed_cuda);

    // start = init_time();
    // int *final_state = cpu_game_of_life(initial_state, n, m, nsteps);
    // float elapsed = get_elapsed_time(start);
    // printf("Finnished GOL in %f ms\n", elapsed);

    // compare(final_state, final_state_cuda, n, m);
    // visualise(VISUAL_ASCII, 1, initial_state, n, m);
    // visualisation2D(game_of_life_grid, n, m);
    // visualise(VISUAL_ASCII, 100, final_state_cuda, n, m);

    // free(final_state_cuda);
    // free(final_state);
    free(initial_state);
    free(game_of_life_grid);

    free(opt);
    return 0;
}
