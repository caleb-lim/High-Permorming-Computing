#include "common.h"

void get_neighbours_1d_x(struct Options *opt, int i, int *n_i) {
    int num_neighbours = 8;
    n_i[0] = i - 1;
    n_i[1] = i - 1;
    n_i[2] = i - 1;
    n_i[3] = i;
    n_i[4] = i + 1;
    n_i[5] = i + 1;
    n_i[6] = i + 1;
    n_i[7] = i;
    if (opt->ineighbourtype == NEIGHBOUR_EXTENDED)
    {
        num_neighbours = 24;
        n_i[8] = i-1;
        n_i[9] = i-2;
        n_i[10] = i-2;
        n_i[11] = i-2;
        n_i[12] = i-2;
        n_i[13] = i-2;
        n_i[14] = i-1;
        n_i[15] = i;
        n_i[16] = i+1;
        n_i[17] = i+2;
        n_i[18] = i+2;
        n_i[19] = i+2;
        n_i[20] = i+2;
        n_i[21] = i+2;
        n_i[22] = i+1;
        n_i[23] = i;
    }
    if (opt->iboundarytype == BOUNDARY_TORAL || opt->iboundarytype == BOUNDARY_TORAL_X_HARD_Y) {
        int period = opt->n;
        for (int i=0;i<num_neighbours;i++) {
            if (n_i[i] < 0) n_i[i] += period;
            if (n_i[i] >= period) n_i[i] -= period;
        }
    }
}

void get_neighbours_1d_y(struct Options *opt, int j, int *n_j) {
    int num_neighbours = 8;
    n_j[0] = j + 1;
    n_j[1] = j;
    n_j[2] = j - 1;
    n_j[3] = j - 1;
    n_j[4] = j - 1;
    n_j[5] = j;
    n_j[6] = j + 1;
    n_j[7] = j + 1;
    if (opt->ineighbourtype == NEIGHBOUR_EXTENDED)
    {
        num_neighbours = 24;
        n_j[8] = j+2;
        n_j[9] = j+2;
        n_j[10] = j+1;
        n_j[11] = j;
        n_j[12] = j-1;
        n_j[13] = j-2;
        n_j[14] = j-2;
        n_j[15] = j-2;
        n_j[16] = j-2;
        n_j[17] = j-2;
        n_j[18] = j-1;
        n_j[19] = j;
        n_j[20] = j+1;
        n_j[21] = j+2;
        n_j[22] = j+2;
        n_j[23] = j+2;
    }
    if (opt->iboundarytype == BOUNDARY_TORAL || opt->iboundarytype == BOUNDARY_TORAL_Y_HARD_X) {
        int period = opt->m;
        for (int i=0;i<num_neighbours;i++) {
            if (n_j[i] < 0) n_j[i] += period;
            if (n_j[i] >= period) n_j[i] -= period;
        }
    }
}

int neighbour_check(struct Options *opt, int n, int m,
    int *n_i, int *n_j, int *current_grid){
    int neighbours = 0;
    int num_neighbours = 8;
    if (opt->ineighbourtype == NEIGHBOUR_EXTENDED) num_neighbours = 24;
    for (int i=0;i<num_neighbours;i++) {
        if (n_i[i] <0 || n_i[i]>= n) continue;
        if (n_j[i] <0 || n_j[i]>= m) continue;
        int index = m*n_i[i]+n_j[i];
        if (current_grid[index] != ALIVE) continue;
        neighbours++;
    }
    return neighbours;
}

void update_grid(struct Options *opt, int index, int neighbours,
    int *current_grid, int *next_grid)
{
    if (opt->iruletype == RULE_STANDARD) {
        if(current_grid[index] == ALIVE && (neighbours == 2 || neighbours == 3)){
            next_grid[index] = ALIVE;
        }else if(current_grid[index] == DEAD && neighbours == 3){
            next_grid[index] = ALIVE;
        }else{
            next_grid[index] = DEAD;
        }
    }
}

void game_of_life(struct Options *opt,
    int *current_grid, int *next_grid, int n, int m){
    int neighbours;
    unsigned long long index;
    int *neighbour_index_i, *neighbour_index_j;
    //default is 9 point stencile
    int num_neighbours = 8;
    // if extended, using 25 point stencile
    if (opt->ineighbourtype == NEIGHBOUR_EXTENDED) num_neighbours = 24;
    neighbour_index_i = (int *)malloc(sizeof(int)*num_neighbours);
    neighbour_index_j = (int *)malloc(sizeof(int)*num_neighbours);
    for(int i = 0; i < n; i++){
        get_neighbours_1d_x(opt, i, neighbour_index_i);
        for(int j = 0; j < m; j++){
            get_neighbours_1d_y(opt, j, neighbour_index_j);
            index = m * i + j;
            neighbours = neighbour_check(opt, n, m,
                neighbour_index_i, neighbour_index_j, current_grid);
            update_grid(opt, index, neighbours, current_grid, next_grid);
        }
    }
    free(neighbour_index_i);
    free(neighbour_index_j);
}

void game_of_life_stats(struct Options *opt, int step, int *current_grid){
    unsigned long long num_in_state[NUMSTATES];
    int m = opt->m, n = opt->n;
    for(int i = 0; i < NUMSTATES; i++) num_in_state[i] = 0;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
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
    unsigned long long nbytes = sizeof(int) * n * m;

    int *grid = (int *) malloc(nbytes);
    int *updated_grid = (int *) malloc(nbytes);
    if(!grid || !updated_grid){
        printf("Error while allocating memory.\n");
        return -1;
    }
    int current_step = 0;
    int *tmp = NULL;
    generate_IC(opt->iictype, grid, n, m);
    struct timeval start;
    start = init_time();
    while(current_step != nsteps){
        struct timeval localstart;
        localstart = init_time();
        game_of_life_stats(opt, current_step, grid);
        visualise(opt->ivisualisetype, current_step, grid, n, m);
        game_of_life(opt, grid, updated_grid, n, m);
        get_elapsed_time(localstart);
        // swap current and updated grid
        tmp = grid;
        grid = updated_grid;
        updated_grid = tmp;
        current_step++;
    }
    printf("Finnished GOL\n");
    get_elapsed_time(start);
    free(grid);
    free(updated_grid);
    free(opt);
    return 0;
}
