#define MAXITER 1000
#define N       1000

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h> 
#include <string.h>

int main(int argc, char *argv[]) {
    int rank, ncpu, x_idx;
    int start, end, cellIdx;
    float *x, *matrix;
    float *finalArray;
    int local_row;
    int which_rank;
    int startIdx;
    float complex z, kappa;
    int i, j, i_1, j_1, k, green, blue;
    double tstart, tend, elapsed;
    FILE *fp;

    /* Basic Initialization */
    MPI_Init(&argc, &argv); 
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &ncpu);

    // if (ncpu != 10) {
    //     printf("Please set 10 processes");
    //     MPI_Finalize();
    //     return 0;
    // }

    if (ncpu > N) {
        if(rank==0){
            printf("Please have less than %d processes\n",ncpu);
        }
        MPI_Finalize();
        return 0;
    }

    /* Each process should calculate the Mandelbrot set for a cyclically assigned range of cells */
    /* Calculates the number of row that should be calculated for each rank*/
    x=(float *) malloc(N*sizeof(float));
    finalArray = (float *) malloc(N*N*sizeof(float));
    
    
    if(rank == 0){
        //Start Time
        tstart = MPI_Wtime();
    }

    //  Cyclic stripe the rows to all the ranks
    for(int i = 0; i < N; i++){
        // Which row in the sub-matrix are we accessing?
        local_row = i / ncpu;
        // Which rank does this row belong to?
        which_rank = i % ncpu;

        if(rank == which_rank){
            startIdx = i * N;
            start = startIdx;
            end = startIdx + N;
            // printf("[Worker %d]: Processing Row %d\tStarting Index: %d\t Ending Index: %d\n",rank,i,start,end);
            int idx = 0;

                for (int cellIdx = start; cellIdx < end; cellIdx++) {

                    i_1 = (cellIdx) % N;
                    j_1 = (cellIdx) / N;

                    z=kappa= (4.0*(i_1-N/2))/N + (4.0*(j_1-N/2))/N * I;

                    k = 1;
                    while ((cabs(z) <= 2) && (k++ < MAXITER)) 
                        z = z * z + kappa;
                    
                    float results = log((float)k) / log((float)MAXITER);
                    x[idx] = results;
                    // printf("idx: %d\t%f\n",cellIdx,results);
                    idx ++; 
                }

            // printf("[INFO: WORKER %d] DATA SENT: ", rank);
            // for (int i = 0; i < N; i++) {
            //     printf("%f ", x[i]);
            // }
            // printf("\n");

            
        if(rank != 0){
            MPI_Send(   x,                      /* message buffer */
                        N,                      /* one data item */
                        MPI_FLOAT,              /* data item is an float */
                        0,                      /* destination process rank */
                        local_row,              /* user chosen message tag */
                        MPI_COMM_WORLD);    
        }
        else{
            for (int idx=0; idx<N; idx++){
                int finalArrayIndex = N*(local_row*ncpu)+idx;
                // printf("%d \n",finalArrayIndex);
                finalArray[finalArrayIndex] = x[idx];    
            }
        }

        // Free the memory for the local row
        // free(x);
        }
    }    

// Barrier to track when calculations are done
MPI_Barrier(MPI_COMM_WORLD);

// Allocate memory for the global matrix on the root process
if (rank == 0) {
    MPI_Status status;
    int iter = (N / ncpu)+1;
    // printf("Iter: %d\n",iter);

    for(int i=0; i<iter; i++){
        for(int j=1; j<ncpu; j++){
            int check = i*ncpu + j;

            // printf("check: %d ",check);

            if(check > N-1){
                printf("Stop!\n");
                break;
            }
            else{
                MPI_Recv(x, N, MPI_FLOAT, j, i, MPI_COMM_WORLD, &status);
                int tag = status.MPI_TAG;
                
                // printf("[INFO: WORKER [%d][%d]] DATA RECIEVED\n", j,i);
                // for (int l = 0; l < N; l++) {
                //     printf("%f ", x[l]);
                // }
                // printf("\n");
            
                for (int idx=0; idx<N; idx++){
                    int finalArrayIndex = N*(j+(tag*ncpu))+idx;
                    // printf("%d \n",finalArrayIndex);
                    finalArray[finalArrayIndex] = x[idx];    
                }
            
            }
        } 
    }
}

    // Print time
    if(rank == 0){
            // Print time
        tend = MPI_Wtime();
        elapsed = tend - tstart;
        printf("Process elapsed time: %.6f seconds\n", elapsed);
    }



// Gather all the local matrices into the global matrix
// MPI_Gather(x, num_rows * N, MPI_FLOAT, global_matrix, num_rows * N, MPI_FLOAT, 0, MPI_COMM_WORLD);

// Free the memory for the local matrix
// free(local_matrix);

// Print the global matrix on the root process
if (rank == 0) {
    printf("[INFO: ROOT %d] Final Matrix:\n", rank);
    printf("Writing mandelbrot.ppm\n");
    fp = fopen("q3-1000x1000.ppm", "w");
    fprintf(fp, "P3\n%4d %4d\n255\n", N, N);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // printf("%f ", finalArray[i * N + j]);

            if (finalArray[i * N + j] < 0.5) {
            green = (int)(2 * finalArray[i * N + j] * 255);
            fprintf(fp, "%3d\n%3d\n%3d\n", 255 - green, green, 0);
            }else {
            blue = (int)(2 * finalArray[i * N + j] * 255 - 255);
            fprintf(fp, "%3d\n%3d\n%3d\n", 0, 255 - blue, blue);
            }
        }
        // printf("\n");
    }
    fclose(fp);      

    free(finalArray);
}


MPI_Finalize();
return 0;
}