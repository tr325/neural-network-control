
const int SIZE = 3;

/* Solves lyapunov eqn for 0 = L(P,A,U,s)   */ 
/* Returns pointer to P                     */
double* Lyap(double A[], double B[]);


/* Returns the smoothed spectral abcissa for W, when given a specific  */ 
/* epsilon to work towards                          */  
double SimpleSSA(double *W[], double eps);


/* Converts the input matrix to Schur form  */
void Schur(double *A[]);


/*  Converts input matrix to Column-major 1-D array for     */
/*  passing to Fortran subroutines                          */
double* FArrayConvert(double *A[], int dimA);

/*  Repopulates array cA with the results of a Fortran  */
/*  subroutine, held in array fA                        */
void CArrayConvert(double fA[], double *cA[], int dimCA);

