/*!

*/
#include "common.h"

/// get some basic timing info
struct timeval init_time(){
    struct timeval curtime;
    gettimeofday(&curtime, NULL);
    return curtime;
}
/// get the elapsed time relative to start, return current wall time
struct timeval get_elapsed_time(struct timeval start){
    struct timeval curtime, delta;
    gettimeofday(&curtime, NULL);
    delta.tv_sec = curtime.tv_sec - start.tv_sec;
    delta.tv_usec = curtime.tv_usec - start.tv_usec;
    double deltas = delta.tv_sec+delta.tv_usec/1e6;
    printf("Elapsed time %f s\n", deltas);
    return curtime;
}
/// get the memory used as reported by system for process
void report_memory_usage(char where[1000]) 
{
    // apple does not provide /proc/self/status file so 
    // just ignore mem usage 
    #ifdef __APPLE__
    return;
    #endif
    struct memusage
    {
        size_t mem[4];
//vm_current, vm_peak, rss_current, rss_peak;
    } usage;

    const char *stat_file = "/proc/self/status";
    FILE *fp;
    char line[10000] = {0};
    char searchstr[4][100] = {"VmSize:", "VmPeak:", "VmRSS:", "VmHWM:"};
    char tmpstr[100];
    int i;
    fp = fopen (stat_file, "r");
    if (! fp) {
        printf("========================================================= \n");
        printf("Couldn't open %s for memory usage reading\n",stat_file);
        printf("========================================================= \n");
        return;
    }
    while (fgets(line, 10000, fp))
    {
	for (i=0;i<4;i++) 
        {
	    if (strstr(line,searchstr[i]))
            {
                sscanf(&line[7],"%zu", &usage.mem[i]);
            }
        }
    }
    fclose(fp);
    printf("========================================================= \n");
    printf("Memory Usage @ %s \n (VM, Peak VM, RSS, Peak RSS) = (%zu, %zu, %zu, %zu) kilobytes \n", 
        where, usage.mem[0], usage.mem[1], usage.mem[2], usage.mem[3]);
    printf("========================================================= \n");
    printf(" \n");
}


#ifdef __APPLE__

static inline void
CPU_ZERO(my_cpu_set_t *cs) { cs->count = 0; }

static inline void
CPU_SET(int num, my_cpu_set_t *cs) { cs->count |= (1 << num); }

static inline int
CPU_ISSET(int num, my_cpu_set_t *cs) { return (cs->count & (1 << num)); }

int sched_getaffinity(pid_t pid, size_t cpu_size, my_cpu_set_t *cpu_set)
{
  int32_t core_count = 0;
  size_t  len = sizeof(core_count);
  int ret = mysysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
  if (ret) {
    printf("error while get core count %d\n", ret);
    return -1;
  }
  cpu_set->count = 0;
  for (int i = 0; i < core_count; i++) {
    cpu_set->count |= (1 << i);
  }
  return 0;
}

//#ifdef __APPLE__
int mysysctlbyname(const char *name, void *oldp, size_t *oldlenp, void *newp,
	     size_t newlen)
{
	int name2oid_oid[2];
	int real_oid[CTL_MAXNAME+2];
	int error;
	size_t oidlen;

	name2oid_oid[0] = 0;	/* This is magic & undocumented! */
	name2oid_oid[1] = 3;

	oidlen = sizeof(real_oid);
	error = sysctl(name2oid_oid, 2, real_oid, &oidlen, (void *)name,
		       strlen(name));
	if (error < 0) 
		return error;
	oidlen /= sizeof (int);
	error = sysctl(real_oid, oidlen, oldp, oldlenp, newp, newlen);
	return (error);
}

#endif 

#ifdef __APPLE__
void cpuset_to_cstr(my_cpu_set_t *mask, char *str)
#else
void cpuset_to_cstr(cpu_set_t *mask, char *str)
#endif
{
  char *ptr = str;
  int i, j, entry_made = 0;
  int32_t core_count;
#ifdef __APPLE__
  size_t len = sizeof(core_count);
  mysysctlbyname(SYSCTL_CORE_COUNT, &core_count, &len, 0, 0);
#else 
  core_count = get_nprocs();
#endif
  for (i = 0; i < core_count; i++) {
    if (CPU_ISSET(i, mask)) {
      int run = 0;
      entry_made = 1;
      for (j = i + 1; j < core_count; j++) {
        if (CPU_ISSET(j, mask)) run++;
        else break;
      }
      if (!run)
        sprintf(ptr, "%d ", i);
      else if (run == 1) {
        sprintf(ptr, "%d,%d ", i, i + 1);
        i++;
      } else {
        sprintf(ptr, "%d-%d ", i, i + run);
        i += run;
      }
      while (*ptr != 0) ptr++;
    }
  }
  ptr -= entry_made;
  ptr = NULL;
}

void report_core_binding()
{
    int rank, size;
#ifdef _MPI
    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    //and this processes' rank is
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else 
    rank=0, size=1;
#endif
    char binding_report[10*1024+256] = {"Core binding report:\n "};
#ifdef __APPLE__
    my_cpu_set_t coremask;
#else
    cpu_set_t coremask;
#endif
    char clbuf[7 * 1024], hnbuf[64];
    memset(clbuf, 0, sizeof(clbuf));
    memset(hnbuf, 0, sizeof(hnbuf));
    (void)gethostname(hnbuf, sizeof(hnbuf));
#ifdef _OPENMP
    #pragma omp parallel shared (binding_report) private(coremask, clbuf) 
#endif
    {
        char result[9 * 1024 + 256],tmpstr[8 * 1024 + 256];
        (void)sched_getaffinity(0, sizeof(coremask), &coremask);
        cpuset_to_cstr(&coremask, clbuf);
        sprintf(tmpstr, "On node %s : ",hnbuf);
        strcat(result,tmpstr);
        sprintf(tmpstr,"MPI Rank %d : ", rank);
        strcat(result,tmpstr);
        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        sprintf(tmpstr," Thread %d : ",thread);
        strcat(result,tmpstr);
        sprintf(tmpstr, " Core affinity = %s\n", clbuf);
        strcat(result,tmpstr);
#ifdef _OPENMP 
        #pragma omp critical 
#endif 
        {
            strcat(binding_report,result);
        }
    }
    printf("========================================================= \n");
    printf("%s",binding_report);
    printf("========================================================= \n");
    printf(" \n");
#ifdef _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
}

void visualise(enum VisualiseType ivisualisetype, int step, int *grid, int n, int m){
    if (ivisualisetype == VISUAL_ASCII) visualise_ascii(step, grid, n, m);
    if (ivisualisetype == VISUAL_PNG) visualise_png(step, grid, n, m);
    else visualise_none(step);
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

void visualise_none(int step){
    printf("Game of Life, Step %d:\n", step);
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
    char where[1000];
    sprintf(where,"Function %s @ %d", __func__, __LINE__);
    report_memory_usage(where);
}



/// report the setup
void report_runtime_state(struct Options *opt) {
    int rank, size;
#ifdef _MPI
    //find out how big the SPMD world is
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    //and this processes' rank is
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
#else 
    rank=0, size=1;
#endif
#ifdef _MPI
    if (rank != 0) {
        MPI_Barrier(MPI_COMM_WORLD);
        return;
    }
#endif
    unsigned long long nbytes = sizeof(int) * opt->n * opt->m;

    printf("========================================================= \n");
    printf("GOL Running with following \n");
    printf("========================================================= \n");
    printf("Requesting grid size of (%d,%d), which requires %f GB \n",
    opt->n, opt->m, nbytes/1024.0/1024.0/1024.0);
    printf("Number of steps %d \n", opt->nsteps);
    printf("IC type %d \n",opt->iictype);
    printf("Visualization type %d \n",opt->ivisualisetype);
    printf("Neighbour rule type %d \n",opt->ineighbourtype);
    printf("Boundary rule type %d \n",opt->iboundarytype);
    printf("========================================================= \n");
    printf(" \n");
#ifdef _MPI
    printf("========================================================= \n");
    printf("Running with MPI \n");
    printf("Comm World contains %d processes \n", size);
    printf("========================================================= \n");
    printf(" \n");
#endif
#ifdef _OPENMP 
    printf("========================================================= \n");
    printf("Running with OpenMP, version %d \n", _OPENMP);
    printf("Maximum number of threads %d \n", omp_get_max_threads());
    printf("========================================================= \n");
    printf(" \n");
#endif
    report_core_binding();
    char where[1000];
    sprintf(where,"Function %s @ %d", __func__, __LINE__);
    report_memory_usage(where);
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
#ifndef USEPNG
    if (opt->ivisualisetype == VISUAL_PNG) {
        printf("PNG visualisation not enabled at compile time, turning off visualisation from now on.\n");
    }
#endif
    report_runtime_state(opt);
}

