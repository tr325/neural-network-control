#ifndef UTILITIES_H
#define UTILITIES_H


#include<fstream>

using namespace std; 


/* Populates a 2D array with the data held in an ascii file */ 
void LoadMat(ifstream &file, double* A[], int matDim);

/* Returns dimension of input matrix from the file  */
int MatSize(ifstream &file);

/*  Performs matrix multiplication of matrices A and B,     
 *  storing the result in C                                 */
void MatMult(double *A[], double *B[], double *C[], int SIZE);


/*  Finds the value of the smoothed spectral abcissa using  */
/*  the Newton-Raphson root-finding technique               */
bool Newton(double *Q[], double *P[], double eps, double sA, double prec, double &s, int SIZE);

/*  Calculates the trace of a matrix                        */
double Trace(double *A[], int dimA);


/*  Finds the maximum value of the diagonals of a matrix    */
double MaxDiag(double *A[], int dimA);


/* Solves lyapunov eqn for 0 = L(X,A,U,s)   */ 
/* Populates matrix X                       */
void Lyap(double *A[], double *X[], bool TRAN, int dimA);


/* Returns the smoothed spectral abcissa for W, when given a specific  */ 
/* epsilon to work towards                          */  
double SimpleSSA(double *W[], double *P[], double *Q[], double *V[], double eps, int SIZE);


/* Converts the input matrix to Schur form  */
void Schur(double *A[], double *VS[], int SIZE);


/*  Converts input matrix to Column-major 1-D array for     */
/*  passing to Fortran subroutines                          */
double* FArrayConvert(double *A[], int dimA);

/*  Repopulates array cA with the results of a Fortran  */
/*  subroutine, held in array fA                        */
void CArrayConvert(double fA[], double *cA[], int dimCA);


/*  Outputs matrix M to file filename for viewing       */
/*  Overloaded for matrices of doubles and ints         */
void OutputMat(char* filename, double *M[], int SIZE);
void OutputMat(char* filename, int *M[], int SIZE);


#endif

