//  Header file for optimisation functions defined in optimise.cpp

#ifndef OPTIMISE_H
#define OPTIMISE_H

/*  Takes a weight matrix as an input, and optimises
 *  the final SIZE/2 columns to minimise the SSA using
 *  gradient descent techniques.                        */
void OptimiseWMat(double *W[], int *B[], double eps, int inhibNum, int SIZE);


/*  Forms the gradient matrix, d(SSA)/dW                */
void FormGradMat(double *gradMat[], double *W[], double *QP[], double *S_Vect[], double *S_Vect_T[], double ssa, int SIZE);


/*  Performs gradient descent on second SIZE/2 columns
 *  of W, based on gradients found in gradMat[]         */
void GradDescent(double *W[], double *gradMat[], double descentRate, int inhibNum, int SIZE);


/*  From an input W matrix, finds the corresponding V matrix. 
 *  w_ij = b_ij*exp(v_ij)   ->    v_ij = ln(w_ij/b_ij)  if bij != 0 */ 
void Reparam(double *W[], int *B[], double *V[], int SIZE);


/*  Recalculates the W matrix after gradient descent manipulation of V  */
void RecalcW(double *W[], int *B[], double *V[], int SIZE);


/*  Forms new synapses if some decay to zero            */
/*  Returns number of synapses which have been reformed */
int ReformSyn(double *V[], int *B[], int inhibNum, int SIZE);


#endif
