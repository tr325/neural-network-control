/*  Header file for generateW.cpp, containing prototypes for matrix generation functions
*   required to create an appropriate W matrix for stabilisation                */

#ifndef GENERATEW_H
#define GENERATEW_H

/*  Main function - populates an empty W matrix with data 
 *  CURRENTLY: populating non-zero elements of W with a const synase strength
 *  ie. excite strengths are constant - inhib strengths initialised to a const value */
void GenerateWMat(double *W[], int *B[], int inhibCols, int SIZE, bool pairwise, bool wholeNet, double kappa);

/*  Initialises the B matrix, using pairwise correlation probabilities laid out in G's paper
 *  Only valid (in this form) for const sparcity (hence only one sparcity passed)
 *  currently skews to antisymmetry (not sure why at the moment)                    */
void initBPairCorr(int *B[], double cMin, double exSparce, double kappa, int exCols, int SIZE, bool wholenet);

/*  Initialises the B matrix WITHOUT pairwise correlations, and with the option
 *  to keep the E->E and E->I connections constant, and only reainitialse the
 *  I-E and I->I sections of the network                                        */
void initBnoCorr(int *B[], double sparce, int exCols, int SIZE, bool wholeNet);

/*  Checks that the input matrix obeys Dale's Law, and populates B matrix
 *  NB: Should only be called once (at start).  The reparameterisation 
 *  should take care of keeping it enforced during the optimization.   */
void EnforceDale(double *W[], int *B[], int inhibNum, int SIZE);

/*  Generates number between 0 and 1 for initialisation of W matrix */
double randx(void);


#endif
