/*
* @Author: Cheng-Hung HSIEH
* @Date:   2017-05-30 20:30:25
* @Last Modified by:   jason
* @Last Modified time: 2017-06-01 08:59:06
* @Describe: MPI Matrix Multiplication
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include "mpi.h"

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"

#define w1 200
#define w2 300
#define w3 500

#define MASTER 0  
#define FROM_MASTER 1 
#define FROM_WORKER 2
//#define CPU 
//#define DEBUG

MPI_Status status;

double  M[w1][w2],   
        N[w2][w3], 
      MxN[w1][w3],
  CPU_MxN[w1][w3];   
   
void initMatrix(void);
void printMatrix(void);
void CPU_matMul(void);
void matMul(int rows);

int main(int argc, char *argv[]){

    int procNum;
    int rankID;      
    int numWorkers;
    int rows;  
    int rowDispatch, rowExtra, rowOffset;
    struct timeval starttime, endtime;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rankID);
    MPI_Comm_size(MPI_COMM_WORLD, &procNum);
    
    numWorkers = procNum - 1;

if (rankID == MASTER) {

#ifdef CPU
        gettimeofday(&starttime, NULL);
        CPU_matMul();
        gettimeofday(&endtime, NULL);
        double executime = 0;
        executime = (endtime.tv_sec - starttime.tv_sec) * 1000.0;
        executime += (endtime.tv_usec - starttime.tv_usec) / 1000.0;
        printf(ANSI_COLOR_CYAN"CPU time: %13lf msec\n"ANSI_COLOR_RESET, executime);
#endif
#ifdef DEBUG
    printf("Number of processes which are working = %d\n", numWorkers);
#endif
    initMatrix();

#ifdef DEBUG
    printMatrix();
#endif

    rowDispatch = w1 / numWorkers;
    rowExtra    = w1 % numWorkers;
    rowOffset   = 0;
#ifdef DEBUG
    printf("Row be dispatched: %d parts and remain row: %d\n\n",rowDispatch, rowExtra);
#endif
    gettimeofday(&starttime, NULL);

    for (int i = 1; i <= numWorkers; i++) { 
      
        rows = (i <= rowExtra) ? rowDispatch + 1 : rowDispatch;

#ifdef DEBUG
        printf("Proc[%d] handling the %d row, rowOffset: %d ",i , rows, rowOffset);
#endif

        MPI_Send(&rowOffset,               1,    MPI_INT, i, FROM_MASTER, MPI_COMM_WORLD);
        MPI_Send(&rows     ,               1,    MPI_INT, i, FROM_MASTER, MPI_COMM_WORLD);
        MPI_Send(&M[rowOffset][0], rows * w2, MPI_DOUBLE, i, FROM_MASTER, MPI_COMM_WORLD);
        MPI_Send(&N              ,   w2 * w3, MPI_DOUBLE, i, FROM_MASTER, MPI_COMM_WORLD);

        rowOffset += rows;

#ifdef DEBUG
        printf("size[rows][0]: %d ", rows * w2);
        printf("size N: %d\n", w2 * w3);
#endif

    }

    for (int i = 1; i <= numWorkers; i++) { 

      MPI_Recv(&rowOffset,               1, MPI_INT,    i, FROM_WORKER, MPI_COMM_WORLD, &status);
      MPI_Recv(&rows,                    1, MPI_INT,    i, FROM_WORKER, MPI_COMM_WORLD, &status);
      MPI_Recv(&MxN[rowOffset][0], rows*w3, MPI_DOUBLE, i, FROM_WORKER, MPI_COMM_WORLD, &status);

    }

#ifdef DEBUG
    printf(ANSI_COLOR_GREEN"MxN[%d][%d]"ANSI_COLOR_RESET,w1,w3);
    for (int i=0; i < w1; i++) {
        printf("\n");
        for (int j=0; j < w3; j++){
            printf("%6.2f   ", MxN[i][j]);
        }
    }
    printf ("\n");
#endif
    gettimeofday(&endtime, NULL);
    double executime = 0;
    executime = (endtime.tv_sec - starttime.tv_sec) * 1000.0;
    executime += (endtime.tv_usec - starttime.tv_usec) / 1000.0;
    printf(ANSI_COLOR_CYAN"MPI time: %13lf msec\n"ANSI_COLOR_RESET, executime);
}  // end of MASTER 


if (rankID > MASTER) {

    MPI_Recv(&rowOffset,       1,    MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
    MPI_Recv(     &rows,       1,    MPI_INT, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
    MPI_Recv(        &M, rows*w2, MPI_DOUBLE, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);
    MPI_Recv(        &N,   w2*w3, MPI_DOUBLE, MASTER, FROM_MASTER, MPI_COMM_WORLD, &status);

    matMul(rows);
      
    MPI_Send(&rowOffset,         1,    MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
    MPI_Send(     &rows,         1,    MPI_INT, MASTER, FROM_WORKER, MPI_COMM_WORLD);
    MPI_Send(      &MxN, rows * w3, MPI_DOUBLE, MASTER, FROM_WORKER, MPI_COMM_WORLD);

    } 

    MPI_Finalize();



    return 0;
} 



void initMatrix(void){

    for (int i = 0; i < w1; i++){
        for (int j = 0; j < w2; j++){
            M[i][j]= rand() % 30;
        }
    }
      
    for (int i = 0; i < w2; i++){
        for (int j=0; j < w3; j++){
            N[i][j]= rand() % 40; 
        }
    }   

}




void matMul(int rows){

    for (int k = 0; k < w3; k++){
        for (int i = 0; i < rows; i++) {
            MxN[i][k] = 0.0;
            for (int j = 0; j < w2; j++){
                MxN[i][k] +=  M[i][j] * N[j][k];
            }
        }
    }

}

void CPU_matMul(void){

    for (int i = 0; i < w1; ++i) {
        for (int j = 0; j < w3; ++j) {
            for (int k = 0; k < w2; ++k) {
                CPU_MxN[i][j] += M[i][k] * N[k][j];
            }
        }
    }


}

