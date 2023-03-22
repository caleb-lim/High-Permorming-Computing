#include "common.h"

#define INCLUDE_CPU_VERSION 
#include "game_of_life.c"
// void game_of_life_step(int *current_grid, int *next_grid, int n, int m){
//     int neighbours;
//     int n_i[8], n_j[8];
//     for(int i = 0; i < n; i++){
//         for(int j = 0; j < m; j++){
//             // count the number of neighbours, clockwise around the current cell.
//             neighbours = 0;
//             n_i[0] = i - 1; n_j[0] = j - 1;
//             n_i[1] = i - 1; n_j[1] = j;
//             n_i[2] = i - 1; n_j[2] = j + 1;
//             n_i[3] = i;     n_j[3] = j + 1;
//             n_i[4] = i + 1; n_j[4] = j + 1;
//             n_i[5] = i + 1; n_j[5] = j;
//             n_i[6] = i + 1; n_j[6] = j - 1;
//             n_i[7] = i;     n_j[7] = j - 1;

//             if(n_i[0] >= 0 && n_j[0] >= 0 && current_grid[n_i[0] * m + n_j[0]] == ALIVE) neighbours++;
//             if(n_i[1] >= 0 && current_grid[n_i[1] * m + n_j[1]] == ALIVE) neighbours++;
//             if(n_i[2] >= 0 && n_j[2] < m && current_grid[n_i[2] * m + n_j[2]] == ALIVE) neighbours++;
//             if(n_j[3] < m && current_grid[n_i[3] * m + n_j[3]] == ALIVE) neighbours++;
//             if(n_i[4] < n && n_j[4] < m && current_grid[n_i[4] * m + n_j[4]] == ALIVE) neighbours++;
//             if(n_i[5] < n && current_grid[n_i[5] * m + n_j[5]] == ALIVE) neighbours++;
//             if(n_i[6] < n && n_j[6] >= 0 && current_grid[n_i[6] * m + n_j[6]] == ALIVE) neighbours++;
//             if(n_j[7] >= 0 && current_grid[n_i[7] * m + n_j[7]] == ALIVE) neighbours++;

//             if(current_grid[i*m + j] == ALIVE && (neighbours == 2 || neighbours == 3)){
//                 next_grid[i*m + j] = ALIVE;
//             } else if(current_grid[i*m + j] == DEAD && neighbours == 3){
//                 next_grid[i*m + j] = ALIVE;
//             }else{
//                 next_grid[i*m + j] = DEAD;
//             }
//         }
//     }
// }


int main(int argc, char **argv)
{
    printf("GPU version");
    return 0;
}
