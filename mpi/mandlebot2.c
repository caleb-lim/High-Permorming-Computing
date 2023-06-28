
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h> 
#include <string.h>
#include <stdbool.h>

#define MAXITER    1000
#define N          8000
#define CHUNKSIZE  100000
#define WORKTAG    1
#define DIETAG     2

void master() {
    //Initialise Variables
    printf("CHUNKSIZE: %d\n",CHUNKSIZE);
    int rank, ncpus, i, index, green, blue, loop;
    long grid_size = 0;
    float *grid, *grid_chunk;
    double  tstart,tend,elapsed;
    FILE   *fp;
    MPI_Comm_size(MPI_COMM_WORLD, &ncpus);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // printf("RANK %d and I am MASTER\n", rank);
    grid = (float *) malloc(N*N*sizeof(float));

    //Start Time
    tstart = MPI_Wtime();
    
    //send first few jobs over
    for(i=1; i<ncpus; i++){
        // printf("SENDING CHUNK TO %d STARTING AT INDEX %d\n", i, grid_size);
        MPI_Send(   &grid_size,              /* message buffer */
                    1,                      /* one data item */
                    MPI_INT,                /* data item is an float */
                    i,                      /* destination process rank */
                    0,                      /* user chosen message tag */
                    MPI_COMM_WORLD);      
        grid_size += CHUNKSIZE;
    }

    // printf("GRID SIZE: %d\n", grid_size);
    int done = 0;           /* counts workers that are done */
    int num_workers = ncpus-1;  /* sent already ncpus */
    
    while (done < num_workers && grid_size < (N*N+10)) {
        int flag = 0;
        MPI_Status status;
        // MPI_Status_init(&status);
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &flag, &status);

        if (flag) {
            int worker_recv = status.MPI_SOURCE;

            // Allocate memory for grid_chunk
            grid_chunk = (float *) malloc((CHUNKSIZE + 1) * sizeof(float));
            MPI_Recv(grid_chunk, CHUNKSIZE + 1, MPI_FLOAT, worker_recv, 0, MPI_COMM_WORLD, &status);

            int starting_index = (int)grid_chunk[0];
            // printf("[INFO: MASTER] DATA RECIEVED FROM WORKER %d: ", worker_recv);
            // for (int l = 0; l < CHUNKSIZE + 1; l++) {
            //     printf("%f ", grid_chunk[l]);
            // }
            // printf("\n");

            for (int k = 1; k <= CHUNKSIZE; k++) {
                grid[starting_index + (k - 1)] = grid_chunk[k];
            }
            

            // free(grid_chunk); 
            // If grid_size < size of matrix then more data will be send to the corresponding worker
            // Else send a DIETAG message to notify the worker to kill the process
            if (grid_size < (N*N)) {
                    // printf("SENDING CHUNK TO WORKER %d STARTING AT INDEX %d\n", worker_recv, grid_size);
                    MPI_Send(   &grid_size,             /* message buffer */
                                1,                      /* one data item */
                                MPI_INT,                /* data item is an float */
                                worker_recv,            /* destination process rank */
                                0,                      /* user chosen message tag */
                                MPI_COMM_WORLD); 
                    grid_size += CHUNKSIZE;

            }
            else {
                    MPI_Send(   
                    &grid_size,             /* message buffer */
                    1,                      /* one data item */
                    MPI_INT,                /* data item is an float */
                    worker_recv,                      /* destination process rank */
                    DIETAG,                 /* user chosen message tag */
                    MPI_COMM_WORLD); 
                    done++;
                    printf("Worker done: %d\n", done);
            }
        }
    }
    
    printf("EXISTING WHILE LOOP\n");

    // Print time
    tend = MPI_Wtime();
    elapsed = tend - tstart;
    printf("Process elapsed time: %.6f seconds\n", elapsed);


    
    // for (int i = 0; i < N * N; i++) {
    //     printf("%f ", grid[i]);
    // if ((i + 1) % N == 0)
    //     printf("\n");
    // }

    /* ----------------------------------------------------------------*/
  
    printf("Writing mandelbrot.ppm\n");
    fp = fopen ("q2-8000x8000.ppm", "w");
    fprintf (fp, "P3\n%4d %4d\n255\n", N, N);
    
    for (loop=0; loop<N*N; loop++) 
	    if (grid[loop]<0.5) {
	        green= (int) (2*grid[loop]*255);
                fprintf (fp, "%3d\n%3d\n%3d\n", 255-green,green,0);
	    } else {
	        blue= (int) (2*grid[loop]*255-255);
                fprintf (fp, "%3d\n%3d\n%3d\n", 0,255-blue,blue);
	    }
    
    fclose(fp);

/* ----------------------------------------------------------------*/
}

void worker() {
    int rank, index, tag;
    int start, end, work_index;
    float *work;
    float complex z, kappa;
    MPI_Status          status;

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    while (true){
        MPI_Recv(   &index, 
                    1, 
                    MPI_INT, 
                    0, 
                    MPI_ANY_TAG, 
                    MPI_COMM_WORLD, 
                    &status);
    tag = status.MPI_TAG;        
    // printf("RANK %d and I am WORKER\n", rank);
    // printf("[INFO: WORKER %d]: TAG: %d\n", rank, tag);
       start = index;
        end = index + CHUNKSIZE;
        // printf("[INFO: WORKER %d] RECIEVE CHUNK FROM MASTER Start: %d End: %d\n",rank,start,end);

        if(tag==0){ 
            work=(float *) malloc((CHUNKSIZE+1)*sizeof(float));
            work[0] = index;
            work_index = 1;
            for (long loop = start; loop < end; loop++) {
                if (loop < N * N) { // Check if the loop index is within the assigned range
                    int i = (loop % N);
                    int j = (loop / N);

                    float complex z, kappa;
                    z = kappa = (4.0 * (i - N / 2)) / N + (4.0 * (j - N / 2)) / N * I;

                    int k = 1;
                    while ((cabs(z) <= 2) && (k++ < MAXITER))
                        z = z * z + kappa;

                    // Process the result as needed
                    float result = log((float)k) / log((float)MAXITER);
                    work[work_index] = result;
                    work_index ++; 
                }
            }

            // printf("[INFO: WORKER %d] SENDING CHUNK TO MASTER Start: %d End: %d\n",rank,start,end);
            // printf("[INFO: WORKER %d] DATA SENT: ", rank);
            // for (int i = 0; i < CHUNKSIZE + 1; i++) {
            //     printf("%f ", work[i]);
            // }
            // printf("\n");

            MPI_Send(   work,                   /* message buffer */
                        CHUNKSIZE+1,            /* one data item */
                        MPI_FLOAT,              /* data item is an float */
                        0,                      /* destination process rank */
                        0,                      /* user chosen message tag */
                        MPI_COMM_WORLD);    
        }
        // else (status.MPI_TAG == DIETAG) 
        else{
            return;
        }
    }
}

int main(int argc, char *argv[]) {
    int rank, ncpu;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpu);

    if (ncpu != 11) {
        if(rank==0)
            printf("Please set 11 processes\n");
        MPI_Finalize();
        return 0;
    }

    // Root node will be initialized as the master and distribute chunk sizes to the worker nodes
    if (!rank) {
        master();
    } else {
        worker();
    }

    MPI_Finalize();
    return 0;
}