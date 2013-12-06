//  Header file for optimisation functions defined in optimise.cpp


/*  Takes a weight matrix as an input, and optimises
 *  the final SIZE/2 columns to minimise the SSA using
 *  gradient descent techniques.                        */
void OptimiseWMat(double *W[], double eps, int SIZE);


/*  Forms the gradient matrix, d(SSA)/dW                */
void FormGradMat(double *gradMat[], double *W[], double *QP[], double ssa, int SIZE);


/*  Performs gradient descent on second SIZE/2 columns
 *  of W, based on gradients found in gradMat[]         */
void GradDescent(double *W[], double *gradMat[], double descentRate, int SIZE);
    
