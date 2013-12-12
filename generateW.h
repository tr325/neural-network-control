/*  Header file for generateW.cpp, containing prototypes for matrix generation functions
*   required to create an appropriate W matrix for stabilisation                */


/*  Main function - populates an empty W matrix with data 
 *  CURRENTLY: populating non-zero elements of W with a const synase strength
 *  ie. excite strengths are constant - inhib strengths initialised to a const value */
void GenerateWMat(double *W[], int *B[], int SIZE);

/*  Checks that the input matrix obeys Dale's Law, and populates B matrix
 *  NB: Should only be called once (at start).  The reparameterisation 
 *  should take care of keeping it enforced during the optimization.   */
void EnforceDale(double *W[], int *B[], int inhibNum, int SIZE);

