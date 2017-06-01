/*
* @Author: Cheng-Hung HSIEH
* @Date:   2017-05-30 20:30:25
* @Last Modified by:   jason
* @Last Modified time: 2017-06-01 08:58:20
* @Describe: MPI Matrix Multiplication
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
//#define DEBUG

int w1, 
    w2, 
    w3;
void initMatrix(double *M, double *N);
void CPU_matMul(double *M, double *N, double *CPU_MxN, int w1, int w2, int w3);

int main(int argc, char *argv[]){


	if(argc != 4 ){
		fprintf(stderr, "Please enter the matrix M[w1][w2] N[w2][w3] CPU_MxN[w1][w3]\n");
	}

	double *M, *N, *CPU_MxN;

	w1 = atoi(argv[1]);
	w2 = atoi(argv[2]);
	w3 = atoi(argv[3]);
	


		  M = (double *) calloc( w1 * w2, sizeof(double));
    	  N = (double *) calloc( w2 * w3, sizeof(double));
    CPU_MxN = (double *) calloc( w1 * w3, sizeof(double));

    initMatrix(M, N);

#ifdef DEBUG
    /*  M[w1][w2] */
    printf(ANSI_COLOR_GREEN"Matrix M[%d][%d]\n"ANSI_COLOR_RESET, w1, w2);
    for (int i = 0; i < w1; i++) {
        for (int j = 0; j < w2; j++) {
            printf("%6.2f", M[i * w2 + j]);
        }
        printf("\n");
    }
    /*  N[w2][w3] */
    printf(ANSI_COLOR_GREEN"Matrix N[%d][%d]\n"ANSI_COLOR_RESET, w2, w3);
    for (int i = 0; i < w2; i++) {
        for (int j = 0; j < w3; j++) {
            printf("%6.2f", N[i * w3 + j]);
        }
        printf("\n");
    }
#endif
	struct timeval starttime, endtime;
	gettimeofday(&starttime, NULL);

    CPU_matMul(M, N, CPU_MxN, w1, w2, w3);

    gettimeofday(&endtime, NULL);
    


#ifdef DEBUG
    printf(ANSI_COLOR_GREEN"CPU_MxN[%d][%d]"ANSI_COLOR_RESET,w1,w3);
    for (int i=0; i < w1; i++) {
        printf("\n");
        for (int j=0; j < w3; j++){
            printf("%6.2f\t", CPU_MxN[i * w3 + j]);
        }
    }
    printf ("\n");
#endif

    double executime = 0;
    executime = (endtime.tv_sec - starttime.tv_sec) * 1000.0;
    executime += (endtime.tv_usec - starttime.tv_usec) / 1000.0;

    printf(ANSI_COLOR_RED"CPU time: %13lf msec\n"ANSI_COLOR_RESET, executime);
	return 0;
}
void CPU_matMul(double *M, double *N, double *CPU_MxN, int w1, int w2, int w3){

    for (int i = 0; i < w1; ++i) {
        for (int j = 0; j < w3; ++j) {
            for (int k = 0; k < w2; ++k) {
                CPU_MxN[i * w3 + j] += M[i * w2 + k] * N[k * w3 + j];
            }
        }
    }


}

void initMatrix(double *M, double *N){

    for (int i = 0; i < w1; i++){
        for (int j = 0; j < w2; j++){
            M[i * w2 + j]= rand() % 30;
        }
    }
      
    for (int i = 0; i < w2; i++){
        for (int j=0; j < w3; j++){
            N[i * w3 + j]= rand() % 40; 
        }
    }   

}