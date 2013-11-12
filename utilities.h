#include<fstream>

using namespace std; 

const int SIZE = 50;


/* Populates a 2D array with the data held in an ascii file */ 
void loadMat(ifstream &file, double* A[], int matDim);


/*  Performs matrix multiplication of matrices A and B,     
 *  storing the result in C                                 */
void MatMult(double A[], double B[], double C[]);


/*  Finds the value of the smoothed spectral abcissa using  */
/*  the Newton-Raphson root-finding technique               */
bool Newton(double *Q[], double *P[], double eps, double sA, double prec, double &s);

/*  Calculates the trace of a matrix                        */
double Trace(double *A[], int dimA);


/*  Finds the maximum value of the diagonals of a matrix    */
double MaxDiag(double *A[], int dimA);


/* Solves lyapunov eqn for 0 = L(P,A,U,s)   */ 
/* Returns pointer to P                     */
double* Lyap(double A[], bool TRAN);


/* Returns the smoothed spectral abcissa for W, when given a specific  */ 
/* epsilon to work towards                          */  
double SimpleSSA(double *W[], double eps);


/* Converts the input matrix to Schur form  */
void Schur(double A[]);


/*  Converts input matrix to Column-major 1-D array for     */
/*  passing to Fortran subroutines                          */
double* FArrayConvert(double *A[], int dimA);

/*  Repopulates array cA with the results of a Fortran  */
/*  subroutine, held in array fA                        */
void CArrayConvert(double fA[], double *cA[], int dimCA);

