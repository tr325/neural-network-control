//  File to test optimization of W matrix, using SSA solving methods in utilities.cpp


#include<iostream>
#include<cmath>
#include"optimise.h"
#include"utilities.h"

using namespace std;


/*  Takes the input weighting matrix, and optimises the second half
 *  to minimise the ssa value.  Optimisation is by gradient descent */
void OptimiseWMat(double *W[], double eps, int SIZE)
{
    double *P[SIZE];
    double *Q[SIZE];
    double *gradMat[SIZE];
    double *PQ[SIZE];
    double *fP;
    double *fQ;
    double *fPQ;
    bool conv;
    double ssa;
    double ssaOld;
    double precision;
    int loopcount;

    for(int i=0; i<SIZE; i++)
    {
        P[i] = new double[SIZE];
        Q[i] = new double[SIZE];
        gradMat[i] = new double[SIZE];
        PQ[i] = new double[SIZE];
    }

    conv = false;
    ssa = 0.0;
    precision = 0.00001;
    loopcount = 0; 
    
    // loops until the gradient descent algorithm has converged
    while(!conv)
    {
        cout << "Optimisation loopcount = " <<loopcount << endl; 
        // solve SSA for current W matrix
        ssaOld = ssa;
        ssa = SimpleSSA(W, P, Q, eps, SIZE);

        // find the gradient matrix, d(SSA)/dW
        fP = FArrayConvert(P, SIZE);
        fQ = FArrayConvert(Q, SIZE);
        fPQ = FArrayConvert(PQ, SIZE);
        MatMult(fP, fQ, fPQ, SIZE);
        
        delete[] fP;
        delete[] fQ;
        CArrayConvert(fPQ, PQ, SIZE);
        delete[] fPQ;
        FormGradMat(gradMat, PQ, SIZE);
        
        // perform gradient descent to optimise W
        GradDescent(W, gradMat, SIZE);

        if(abs(ssa - ssaOld) < precision)
        {
            conv = true;
        }
        
        loopcount++;
    }

    for(int i=0; i<SIZE; i++)
    {
        delete[] P[i];
        delete[] Q[i];
        delete[] gradMat[i];
        delete[] PQ[i];
    }

    return;
}


/*  Forms the gradient matrix           */
void FormGradMat(double *gradMat[], double *PQ[], int SIZE)
{
    double trace;

    trace = Trace(PQ, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            gradMat[i][j] = PQ[i][j]/trace;
        }
    }

    return;
}


/*  Performs gradient descent on W, based on gradients in gradMat  */
void GradDescent(double *W[], double *gradMat[], int SIZE)
{
    int descentRate;

    // Determines how far down the gradient each iteration moves
    descentRate = 1;

    for(int i=(SIZE/2); i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            if(i != j)
            {
                W[i][j] -= descentRate*gradMat[i][j];
            }
        }
    }
    
    return;
}












