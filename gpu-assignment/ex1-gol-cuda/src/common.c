/*!

*/
#include "common.h"
void compare(int *cpu, int *gpu, int n, int m){
    int match = 1;
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            if (cpu[i*n + j] != gpu[i*n + j]) {
                match = 0;
                int error = i*m + j;
                printf("Error at %d", error);
            }
        }
    }

    if(match == 1){
         printf("Both matrixes are matching\n");
    }
}

void convert2D(int **grid, int *array, int n, int m){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            grid[i][j] = array[i*m + j];
        }
    }
}


void visualisation2D(int **grid, int n, int m){
    printf("Game of Life\n");
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            char cell = ' ';
            if (grid[i][j] == ALIVE) cell = '*';
            printf(" %c ", cell);
        }
        printf("\n");
    }
}

void visualise(enum VisualiseType ivisualisetype, int step, int *grid, int n, int m){
    if (ivisualisetype == VISUAL_ASCII) visualise_ascii(step, grid, n, m);
    if (ivisualisetype == VISUAL_PNG) visualise_png(step, grid, n, m);
    else visualise_none();
}

/// ascii visualisation
void visualise_ascii(int step, int *grid, int n, int m){
    printf("Game of Life\n");
    printf("Step %d:\n", step);
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {
            char cell = ' ';
            if (grid[i*m + j] == ALIVE) cell = '*';
            printf(" %c ", cell);
        }
        printf("\n");
    }
}

void visualise_png(int step, int *grid, int n, int m){
#ifdef USEPNG
    char pngname[2000];
    sprintf(pngname,"GOL.grid-%d-by-%d.step-%d.png",n,m,step);
    bitmap_t gol;
    gol.width = n;
    gol.height = m;
    gol.pixels = calloc (n*m, sizeof (pixel_t));
    if (! gol.pixels) {
        exit(9);
    }
    for(int i = 0; i < n; i++)
    {
        for(int j = 0; j < m; j++)
        {

            pixel_t * pixel = pixel_at (&gol, i, j);
            int state = grid[i*m+j];
            if (state == ALIVE) {
                pixel->red = (uint8_t)0;
                pixel->green = (uint8_t)255;
                pixel->blue = (uint8_t)0;
            }
            else if (state == DEAD) {
                pixel->red = (uint8_t)0;
                pixel->green = (uint8_t)0;
                pixel->blue = (uint8_t)0;
            }
            else if (state == BORN) {
                pixel->red = (uint8_t)0;
                pixel->green = (uint8_t)255;
                pixel->blue = (uint8_t)255;
            }
            else if (state == DYING) {
                pixel->red = (uint8_t)255;
                pixel->green = (uint8_t)0;
                pixel->blue = (uint8_t)0;
            }
        }
    }

    if (save_png_to_file (&gol, pngname)) {
        fprintf (stderr, "Error writing png file %s\n", pngname);
        exit(9);
    }
    free (gol.pixels);
#endif
}

void visualise_none(){
}

/// generate random IC
void generate_rand_IC(int *grid, int n, int m){
    for(int i = 0; i < n; i++){
        for(int j = 0; j < m; j++){
            grid[i*m + j] = (rand() % 100 < 40) ? DEAD : ALIVE;
        }
    }
}
/// generate some ICs
void generate_IC(enum ICType ic_choice, int *grid, int n, int m){
    if (ic_choice == IC_RAND) generate_rand_IC(grid, n, m);
}

/// get some basic timing info
struct timeval init_time(){
    struct timeval curtime;
    gettimeofday(&curtime, NULL);
    return curtime;
}
/// get the elapsed time relative to start, return current wall time
float get_elapsed_time(struct timeval start){
    struct timeval curtime, delta;
    gettimeofday(&curtime, NULL);
    timersub(&curtime, &start, &delta);
    float elapsed = delta.tv_sec * 1000.0f + delta.tv_usec / 1000.0f;
    return elapsed;
}


/// UI
void getinput(int argc, char **argv, struct Options *opt){
  if(argc < 3){
      printf("Usage: %s <grid height> <grid width> [<nsteps> <IC type> <Visualisation type> <Rule type> <Neighbour type> <Boundary type> <stats filename> ]\n", argv[0]);
      exit(0);
  }
  // grid size
  char statsfilename[2000] = "GOL-stats.txt";
  opt->n = atoi(argv[1]), opt->m = atoi(argv[2]);
  opt->nsteps = -1;
  if (argc >= 4)
      opt->nsteps = atoi(argv[3]);
  if (argc >= 5)
      opt->iictype = atoi(argv[4]);
  if (argc >= 6)
      opt->ivisualisetype = atoi(argv[5]);
  if (argc >= 7)
      opt->iruletype = atoi(argv[6]);
  if (argc >= 8)
      opt->ineighbourtype = atoi(argv[7]);
  if (argc >= 9)
      opt->iboundarytype = atoi(argv[8]);
  if (argc >= 10)
    strcpy(statsfilename, argv[9]);
  if (opt->n <= 0 || opt->m <= 0) {
      printf("Invalid grid size.\n");
      exit(1);
  }
  strcpy(opt->statsfile, statsfilename);
  unsigned long long nbytes = sizeof(int) * opt->n * opt->m;
  printf("Requesting grid size of (%d,%d), which requires %f GB \n",
  opt->n, opt->m, nbytes/1024.0/1024.0/1024.0);
#ifndef USEPNG
  if (opt->ivisualisetype == VISUAL_PNG) {
      printf("PNG visualisation not enabled at compile time, turning off visualisation from now on.\n");
  }
#endif
}
