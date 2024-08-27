#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <mpi.h>


int main(int argc, char *argv[]){

    int N = atoi(argv[1]);
    MPI_Init(NULL, NULL);
    double start_time = MPI_Wtime();
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    
    if (N < 2) {
        exit(1);
    }
    
    char filename[50];
    sprintf(filename, "%d.txt", N);
    
    int num_perProcess = ceil((N+1)/world_size);
    int start = world_rank * num_perProcess + 2;
    int end = start + num_perProcess-1;
    if (world_rank == world_size - 1) {
        end = N; // Ensure the last process covers the range up to 100
    }
    
    if(world_rank!=0){
        int* arr = calloc(N, sizeof(int));
        for (int j = start; j <= end; j ++) {
            arr[j-2]=j;
            for (int i = 2;i <= floor((N + 1) / 2);i++) {
                if(j%i==0&&j!=i){
                    arr[j-2]=0;
                }
            }
        }
        MPI_Send(arr, N, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    if(world_rank==0&&world_size==1){
        int* result = calloc(N, sizeof(int));
        for (int j = 2; j <= N; j ++) {
            result[j-2]=j;
            for (int i = 2;i <= floor((N + 1) / 2);i++) {
                if(j%i==0&&j!=i){
                    result[j-2]=0;
                }
            }
        }
        FILE *file = fopen(filename, "w");
        
        for(int i=0;i<N;i++){
            if(result[i]!=0){
            fprintf(file, "%d ", result[i]);
            }
        }
        double end_time = MPI_Wtime();
        fprintf(file,"\ntime = %f seconds\n",end_time-start_time);
        fclose(file);
    }

    if(world_rank==0&&world_size!=1){
        int* result = calloc(N, sizeof(int));
        int** prime = (int**) malloc(world_size * sizeof(int*));
        if (prime == NULL) {
            exit(1);
        }

        for (int i = 0; i < world_size; i++) {
            prime[i] = (int*) calloc(N, sizeof(int));
            if (prime[i] == NULL) {
                for (int j = 0; j < i; j++) {
                    free(prime[j]);
                }
                free(prime);
                exit(1);
            }
        }
        for (int j = 2; j <= end; j ++) {
            result[j-2]=j;
            for (int i = 2;i <=floor((N + 1) / 2) ;i++) {
                if(j%i==0&&j!=i){
                    result[j-2]=0;
                }
            }
        }
    
        for(int src = 1; src < world_size; src++){
            MPI_Recv(prime[src-1], N, MPI_INT, src, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    
    
        for(int i=0;i<world_size-1;i++){
            for(int j=0;j<N;j++){
                if(prime[i][j] != 0){
                    result[j]=prime[i][j];
                }
            }    
        }
    
        FILE *file = fopen(filename, "w");
        
        for(int i=0;i<N;i++){
            if(result[i]!=0){
            fprintf(file, "%d ", result[i]);
            }
        }
        double end_time = MPI_Wtime();
        fprintf(file,"\ntime = %f seconds\n",end_time-start_time);
        fclose(file);
    }

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Finalize();
    return 0;
}