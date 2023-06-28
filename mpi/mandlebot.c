#define MAXITER 1000
#define N	8000

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <mpi.h> 
#include <string.h>


int main(int argc, char *argv[]) {
    //Initialise MPI variables 
    int     rank,ncpu,numCells,numIter;
    float   *x;
    float   complex z, kappa;
    int     i,j,k,green,blue;
    int start,end,loop;
    double  tstart,tend,elapsed;
    FILE    *fp;

    /* Basic Initialisation */
    MPI_Init(&argc,&argv); 
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&ncpu);

    if(ncpu!=10){
        printf("Please set 10 processes");
        MPI_Finalize();
        return 0;
    }

    // if((N*N)%10!=0){
    //     printf("Can't split it to 10 even segments");
    //     MPI_Finalize();
    //     return 0;
    // }

    /* Each process should split the array into segements of N/10 rows */
    /* Split the grid into 10 equal segments */
    numCells = (N*N)/ncpu;
    

    // Allocated variable x array for each process.
    // However process 9 should take the remaining cells if N*N is not divisible by 10
    if(rank!=9){
        x=(float *) malloc(numCells*sizeof(float));
        // numIter = numCells;
    }
    else{
        x=(float *) malloc((N*N - (numCells*9))*sizeof(float));
        // numIter = N*N - (numCells*9);
    }

    // Set all values in x to 0
    memset(x, 0, numCells * sizeof(float));

    // Now we need to find the start and end index for each process for the for loop
    start = rank * numCells;
    end = start + numCells;
    // printf("RANK %d and I am calculating %d cells. Start: %d End %d\n",rank,(end-start),start,end);

    // Collect start time in rank 0
    if (rank == 0) {
        tstart = MPI_Wtime();
    }

    // We can calculate the mandelbrot set for each process
    for(loop=0; loop<numCells; loop++){
        // Check if the loop index is within the assigned range
        if (loop < N * N) { 
            i=(start+loop)%N;
            j=(start+loop)/N;

            z=kappa= (4.0*(i-N/2))/N + (4.0*(j-N/2))/N * I;

            k=1;
            while ((cabs(z)<=2) && (k++<MAXITER)) 
                z= z*z + kappa;
            
            x[loop]= log((float)k) / log((float)MAXITER);
            // x[loop] = 1;
            // printf("Rank: %d Cell: %d Ans: %f\n",rank,loop,x[loop]);
        }
    }

    // Time it takes to complete one process
    double processTime = MPI_Wtime() - tstart;
    printf("Process elapsed time for RANK %d: %.6f seconds\n",rank,processTime);

    // Add an MPI barrier to synchronize all processes
    // printf("Finished Process Rank %d\n",rank);
    MPI_Barrier(MPI_COMM_WORLD);

    // Use MPI Gather get the processes answer the a final array
    float *finalArray;
    if(rank==0){
        finalArray = (float *) malloc(N*N*sizeof(float));
    }
    MPI_Gather(x,numCells,MPI_FLOAT,finalArray,numCells,MPI_FLOAT,0,MPI_COMM_WORLD);

    //Time it takes to complete all processes
    if(rank == 0){    
        tend = MPI_Wtime();
        elapsed = tend - tstart;
        printf("Process elapsed time: %.6f seconds\n", elapsed);
        fflush(stdout);
    }

    // Print out the final array
    if(rank==0){
        // printf("Final Array\n");
        // for(i=0;i<N*N;i++){
        //     printf("%f ",finalArray[i]);
        //     if(i%10==(N-1)){
        //         printf("\n");
        //     }
        // }
        int print;

        printf("Writing mandelbrot.ppm\n");
        fp = fopen ("q1-8000x8000.ppm", "w");
        fprintf (fp, "P3\n%4d %4d\n255\n", N, N);
        for (print=0; print<N*N; print++) {
            if (finalArray[print]<0.5) {
                green= (int) (2*finalArray[print]*255);
                    fprintf (fp, "%3d\n%3d\n%3d\n", 255-green,green,0);
            } else {
                blue= (int) (2*finalArray[print]*255-255);
                    fprintf (fp, "%3d\n%3d\n%3d\n", 0,255-blue,blue);
            }
        }
        fclose(fp);
    }

    free(x);
    MPI_Finalize();
    return 0;
}