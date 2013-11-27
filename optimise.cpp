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
    double *QP[SIZE];
    double *fP;
    double *fQ;
    double *fQP;
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
        QP[i] = new double[SIZE];
    }

    conv = false;
    ssa = 0.0;
    precision = 0.00001;
    loopcount = 0; 
    
    ofstream resFile;
    resFile.open("SsaOptimise.ascii");
    
    // loops until the gradient descent algorithm has converged
    while(!conv)
    {
        //cout << "Optimisation loopcount = " <<loopcount << endl; 
        // solve SSA for current W matrix
        ssaOld = ssa;
        ssa = SimpleSSA(W, P, Q, eps, SIZE);

        // find the gradient matrix, d(SSA)/dW
       /* fP = FArrayConvert(P, SIZE);
        fQ = FArrayConvert(Q, SIZE);
        fQP = FArrayConvert(QP, SIZE);
        MatMult(fQ, fP, fQP, SIZE);
            
        delete[] fP;
        delete[] fQ;
        CArrayConvert(fQP, QP, SIZE);
        delete[] fQP;*/
        FormGradMat(gradMat, P, Q, W, QP, ssa, SIZE);
        
        ofstream grFile;
        grFile.open("TESTGradMat.ascii");
        for(int i=0; i<SIZE; i++)
        {
            for(int j=0; j<SIZE; j++)
            {
                grFile << gradMat[i][j] << "  ";
            }
            grFile << endl; 
        }                
        
        // perform gradient descent to optimise W
        GradDescent(W, gradMat, SIZE);

        ofstream opfile; 
        opfile.open("optW.ascii");
        for(int i=0; i<SIZE; i++)
        {
            for(int j=0; j<SIZE; j++)
            {
                opfile << W[i][j] << "  ";
            }
            opfile <<endl; 
        }
        
        if(abs(ssa - ssaOld) < precision)
        {
            cout << "converged" <<endl;
            conv = true;
        }
        
        resFile << ssa << endl; 
        
        cout << "SSA on loop " << loopcount << ": " <<ssa <<endl; 
        loopcount++;
    }

    for(int i=0; i<SIZE; i++)
    {
        delete[] P[i];
        delete[] Q[i];
        delete[] gradMat[i];
        delete[] QP[i];
    }

    return;
}


/*  Forms the gradient matrix           */
void FormGradMat(double *gradMat[], double *P[], double *Q[], double *A[], double *QP[], double ssa, int SIZE)
{
    //double *A[SIZE];
   // double *I[SIZE];
    double *fQ;
    double *fP;
    double *fQP;
    double *fA;
    double tr;
 /*  
    for(int i=0; i<SIZE; i++)
    {
        A[i] = new double[SIZE];
        I[i] = new double[SIZE];
    }
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            if(i == j)
            {
                I[i][j] = 1;
            }
            else
            {
                I[i][j] = 0;
            }
            A[i][j] = W[i][j] - ssa*I[i][j];
        }
    }
    */
   // cout << "array exists here " << W[10][10] <<endl; 
    
    fA = FArrayConvert(A, SIZE);
    fQ = Lyap(fA, true, SIZE);
    CArrayConvert(fQ, Q, SIZE);
    
    fP = Lyap(fA, false, SIZE);
    fQP = FArrayConvert(QP, SIZE);
    MatMult(fQ, fP, fQP, SIZE);
    CArrayConvert(fQP, QP, SIZE);
    
    tr = Trace(QP, SIZE);
    cout << "trace(QP) = " <<tr <<endl;
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            gradMat[i][j] = QP[i][j]/tr;
        }
    }
    
    cout <<"exit of FormGradMat()" <<endl; 
   // cout << "does array exist here?" << W[10][10] <<endl; 

    //return;
}


/*  Performs gradient descent on W, based on gradients in gradMat  */
void GradDescent(double *W[], double *gradMat[], int SIZE)
{
    double descentRate;

    // Determines how far down the gradient each iteration moves
    descentRate = 0.02;

    for(int i=(SIZE/2); i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            if(i != j)
            {
                //cout << "Wij updated from " <<W[i][j];
                W[i][j] += -descentRate*gradMat[i][j];
                //cout << " to " <<W[i][j] <<endl;
            }
        }
    }
    
    return;
}












