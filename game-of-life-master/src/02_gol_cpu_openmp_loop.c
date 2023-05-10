#include "common.h"
#include "omp.h"

void game_of_life(struct Options *opt, int *current_grid, int *next_grid, int n, int m){
    int neighbours;
    int n_i[8], n_j[8];

    #pragma omp parallel for private(neighbours, n_i, n_j) collapse(2) schedule(static) shared(current_grid, next_grid)
    for(int i = 0; i < m; i++){
        for(int j = 0; j < n; j++){
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

            //printf("%d of %d\n", omp_get_thread_num(), omp_get_max_threads());  
            

            #pragma omp simd reduction(+:neighbours)       
            { 
            if(n_i[0] >= 0 && n_j[0] >= 0 && current_grid[n_i[0] * m + n_j[0]] == ALIVE) neighbours++;
            if(n_i[1] >= 0 && current_grid[n_i[1] * m + n_j[1]] == ALIVE) neighbours++;
            if(n_i[2] >= 0 && n_j[2] < m && current_grid[n_i[2] * m + n_j[2]] == ALIVE) neighbours++;
            if(n_j[3] < m && current_grid[n_i[3] * m + n_j[3]] == ALIVE) neighbours++;
            if(n_i[4] < n && n_j[4] < m && current_grid[n_i[4] * m + n_j[4]] == ALIVE) neighbours++;
            if(n_i[5] < n && current_grid[n_i[5] * m + n_j[5]] == ALIVE) neighbours++;
            if(n_i[6] < n && n_j[6] >= 0 && current_grid[n_i[6] * m + n_j[6]] == ALIVE) neighbours++;
            if(n_j[7] >= 0 && current_grid[n_i[7] * m + n_j[7]] == ALIVE) neighbours++;
            }


            // for(int k=0; k<8; k++){
            //     if(n_i[k] >= 0 && n_i[k] < n && n_j[k] >= 0 && n_j[k] < m && current_grid[n_i[k] * m + n_j[k]] == ALIVE){
            //         neighbours++;
            //     }
            // }
            
            if(current_grid[i*m + j] == ALIVE && (neighbours == 2 || neighbours == 3)){
                next_grid[i*m + j] = ALIVE;
            } else if(current_grid[i*m + j] == DEAD && neighbours == 3){
                next_grid[i*m + j] = ALIVE;
            }else{
                next_grid[i*m + j] = DEAD;
            }

        }
    }

    // int neighbor_offsets[] = {-m-1,-m,-m+1,-1,1,m-1,m,m+1};
    // for(int i = 0; i < m; i++){
    //     for(int j = 0; j < n; j++){
    //         // count the number of neighbors
    //         int neighbors = 0;
    //         for(int k = 0; k < 8; k++){
    //             int neighbor_index = i*j + neighbor_offsets[k];
    //             if(neighbor_index >= 0 && neighbor_index < m*n && current_grid[neighbor_index] == ALIVE){
    //                 neighbors++;
    //             }
    //         }

    //         // update the next grid
    //         if(current_grid[i*n + j] == ALIVE && (neighbors == 2 || neighbors == 3)){
    //             next_grid[i*n + j] = ALIVE;
    //         } else if(current_grid[i*n + j] == DEAD && neighbors == 3){
    //             next_grid[i*n + j] = ALIVE;
    //         } else {
    //             next_grid[i*n + j] = DEAD;
    //         }
    //     }
    // }

}

void game_of_life_stats(struct Options *opt, int step, int *current_grid){
    unsigned long long num_in_state[NUMSTATES];
    int m = opt->m, n = opt->n;
    for(int i = 0; i < NUMSTATES; i++) num_in_state[i] = 0;
    for(int j = 0; j < m; j++){
        for(int i = 0; i < n; i++){
            num_in_state[current_grid[i*m + j]]++;
        }
    }
    double frac, ntot = opt->m*opt->n;
    FILE *fptr;
    if (step == 0) {
        fptr = fopen(opt->statsfile, "w");
    }
    else {
        fptr = fopen(opt->statsfile, "a");
    }
    fprintf(fptr, "step %d : ", step);
    for(int i = 0; i < NUMSTATES; i++) {
        frac = (double)num_in_state[i]/ntot;
        fprintf(fptr, "Frac in state %d = %f,\t", i, frac);
    }
    fprintf(fptr, " \n");
    fclose(fptr);
}



int main(int argc, char **argv)
{
    struct Options *opt = (struct Options *) malloc(sizeof(struct Options));
    getinput(argc, argv, opt);
    int n = opt->n, m = opt->m, nsteps = opt->nsteps;
    int *grid = (int *) malloc(sizeof(int) * n * m);
    int *updated_grid = (int *) malloc(sizeof(int) * n * m);
    if(!grid || !updated_grid){
        printf("Error while allocating memory.\n");
        return -1;
    }
    int current_step = 0;
    int *tmp = NULL;
    generate_IC(opt->iictype, grid, n, m);
    struct timeval start, steptime;
    start = init_time();
    while(current_step != nsteps){
        steptime = init_time();
        visualise(opt->ivisualisetype, current_step, grid, n, m);
        game_of_life_stats(opt, current_step, grid);
        game_of_life(opt, grid, updated_grid, n, m);
        // swap current and updated grid
        tmp = grid;
        grid = updated_grid;
        updated_grid = tmp;
        current_step++;
        get_elapsed_time(steptime);
    }
    printf("Finnished GOL\n");
    get_elapsed_time(start);

    printf("Running with OpenMP and %d threads\n", omp_get_max_threads());

    free(grid);
    free(updated_grid);
    free(opt);
    return 0;
}
