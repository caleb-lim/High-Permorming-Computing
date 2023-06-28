#include <iostream>
#include <cstring>
#include <complex>
#include <cmath>
#include <mpi.h>

#define MAXITER 1000
#define N	10

int main(int argc, char *argv[]){

    // Declate variables for timing
    double t_start;
    double t_end;
    double t_total;
    int i,j,k,green,blue;


    // Unique rank for this process
    int rank;
    
    // Total number of ranks
    int size;

    // Initializes the MPI execution environment
    MPI_Init(&argc, &argv);
    
    // Get the rank
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the total number ranks in this communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // Calulate the number of rows based on the number of ranks
    int num_rows = N / size; 

    /*
     * Distribute Work to Ranks:
     * Rank 0 needs to send the appropriate rows to each process
     * before they are able to proceed
     */
    // Declare our problem matrices
    // This work is duplicated just for code simplicity
    float *matrix;
    if(rank == 0){
        // Only rank 0 needs space for the total solution 
        matrix = new float[N * N];
    }

    // // Declare our sub-matrix for each process
    // float *sub_matrix = new float[N * num_rows];

    // Cyclic stripe the rows to all the ranks
    if(size == 1){
        // If there is only one rank, we can just copy the matrix
        memcpy(sub_matrix, matrix, N * N * sizeof(float));
    }
    // else{
    //     // If there are multiple ranks, we need to stripe the rows
    //     for(int i = 0; i < num_rows; i++){
    //         // Copy the row into our sub-matrix
    //         memcpy(&sub_matrix[i * N], &matrix[i * size * N], N * sizeof(float));
    //     }
    // }

    /*
     * One rank normalizes the pivot row, then sends it to all
     * later ranks for elimination
     */

    // Get start time
    if(rank == 0){
        t_start = MPI_Wtime();
    }
    
    // Local variables for code clarity
    int local_row;
    int which_rank;
    int startIdx;
    int loop;
    int scale;
    int i_1,j_1;
    int start, end;

    std::complex<float> z, kappa;

    // Allocate space for a single row to be sent to this rank
    float *row = new float[N];

    // Iterate over all rows
    for(int i = 0; i < N; i++){
            
        // Which row in the sub-matrix are we accessing?
        local_row = i / size;
        // Which rank does this row belong to?
        which_rank = i % size;
        
        if(rank == which_rank){
            startIdx = i * N;
            start = startIdx;
            end = startIdx + N;
            std::cout << "[Worker "<<rank<<"]: Processing Row #" << i << "\tStarting Index:"<< start << "\tEnding Index:"<< end <<std::endl;
          
            for(loop=start; loop<end; loop++){
                i_1 = (loop) % N;
                j_1 = (loop) / N;


                std::complex<float> z, kappa;
                z = kappa = (4.0f * ( i_1 - N / 2)) / N + (4.0f * (j_1  - N / 2)) / N * std::complex<float>(0, 1);

                int k = 1;
                while ((std::abs(z) <= 2) && (k++ < MAXITER))
                    z = z * z + kappa;

                float result = log((float)k) / log((float)MAXITER);
                sub_matrix[loop] = result;
                // std::cout << "Rank " << rank << std::endl;
                // std::cout << "Cell " << (loop) << ": " << result << " " << std::endl;
                // x[loop] = 1;
                // printf("Rank: %d Cell: %d Ans: %f\n",rank,loop,x[loop]);
            
            }
            // Copy the row into our send buffer
            memcpy(row, &sub_matrix[startIdx], N * sizeof(float));

            //Print the variable row to see if it is correct
            for(int i = 0; i < N; i++){
                std::cout << row[i] << " ";
            }
            std::cout << std::endl;
        }

        // Send the row to all other ranks
        MPI_Bcast(row, N, MPI_FLOAT, which_rank, MPI_COMM_WORLD);
       
    }

    // Barrier to track when calculations are done
    MPI_Barrier(MPI_COMM_WORLD);

    // Stop the time before the gather phase
    if(rank == 0){
        t_end = MPI_Wtime();
        t_total = t_end - t_start;
    }

    /*
     * Collect all Sub-Matrices
     * All sub-matrices are gathered using the gather function
     */
    // Gather all the sub-matrices into the matrix
    if(rank == 0){

    }
    MPI_Gather(sub_matrix, N , MPI_FLOAT, matrix, N , MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    // Print the matrix
    if(rank == 0){
        for(int i = 0; i < N; i++){
            for(int j = 0; j < N; j++){
                std::cout << matrix[i * N + j] << " ";
            }
            std::cout << std::endl;
        }
    }


    MPI_Finalize();

    // Print the output and the time
    if(rank == 0){
        //print_matrix(matrix, N);
        std::cout << t_total << " Seconds" << std::endl;
    }

    // // Free heap-allocated memory
    // if(rank == 0){
    //     delete[] matrix;
    // }
    // delete[] sub_matrix;
    // delete[] row;

    return 0;
}
