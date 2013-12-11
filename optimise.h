//  Header file for optimisation functions defined in optimise.cpp


/*  Takes a weight matrix as an input, and optimises
 *  the final SIZE/2 columns to minimise the SSA using
 *  gradient descent techniques.                        */
void OptimiseWMat(double *W[], double eps, int SIZE);


/*  Forms the gradient matrix, d(SSA)/dW                */
void FormGradMat(double *gradMat[], double *A[], double *W[], double *QP[], double ssa, int SIZE);


/*  Performs gradient descent on second SIZE/2 columns
 *  of W, based on gradients found in gradMat[]         */
void GradDescent(double *W[], double *gradMat[], double descentRate, int inhibNum, int SIZE);
    
/*  Checks that the input matrix obeys Dale's Law, and populates B matrix
 *  NB: Should only be called once (at start).  The reparameterisation 
 *  should take care of keeping it enforced during the optimization.   */
void EnforceDale(double *W[], int *B[], int inhibNum, int SIZE);


/*  From an input W matrix, finds the corresponding V matrix. 
 *  w_ij = b_ij*exp(v_ij)   ->    v_ij = ln(w_ij/b_ij)  if bij != 0 */ 
void Reparam(double *W[], int *B[], double *V[], int SIZE);


/*  Recalculates the W matrix after gradient descent manipulation of V  */
void RecalcW(double *W[], int *B[], double *V[], int SIZE);
