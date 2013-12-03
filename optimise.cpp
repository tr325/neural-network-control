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
    double *A[SIZE];
    double *I[SIZE];
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
    double descentRate;
    int loopcount;

    for(int i=0; i<SIZE; i++)
    {
        A[i] = new double[SIZE];
        I[i] = new double[SIZE];
        P[i] = new double[SIZE];
        Q[i] = new double[SIZE];
        gradMat[i] = new double[SIZE];
        QP[i] = new double[SIZE];
        for(int j=0; j<SIZE; j++)
        {
            if(i==j)
            {
                I[i][j] = 1;
            }
            else
            {
                I[i][j] = 0;
            }
        }
    }

    conv = false;
    ssa = 0.0;
    ssaOld = 10.0;  //initialised to greater than ssa for descentRate manipulation
    precision = 0.00001;
    loopcount = 0;
    descentRate = 0.5;
    
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
        for(int i=0; i<SIZE; i++)
        {
            for(int j=0; j<SIZE; j++)
            {
                A[i][j] = W[i][j] - ssa*I[i][j];
            }
        }

        fP = FArrayConvert(P, SIZE);
        fQ = FArrayConvert(Q, SIZE);
        fQP = FArrayConvert(QP, SIZE);
        MatMult(fQ, fP, fQP, SIZE);
        CArrayConvert(fQP, QP, SIZE);
        FormGradMat(gradMat, P, Q, A, QP, ssa, SIZE);
        
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
            cout << "Converged" <<endl;
            conv = true;
        }
        /*
        if(loopcount == 100)
        {
            cout <<"That's enough i think..." <<endl;
            conv = true;
        }
        */
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
    double *V[SIZE];
    double *Vt[SIZE];
    double *G[SIZE];
    double *foo[SIZE];
    double *fQ;
    double *fP;
    double *fQP;
    double *fA;
    double *fV;
    double *fVt;
    double *fG;
    double *fFoo;
    double *fGrad;
    double tr;
    
    for(int i=0; i<SIZE; i++)
    {
        V[i] = new double[SIZE];
        Vt[i] = new double[SIZE];
        G[i] = new double[SIZE];
        foo[i] = new double[SIZE];
    }
    
    //  Calculate the matrix of schur vectors, V
    fV = FArrayConvert(V, SIZE);
    fA = FArrayConvert(A, SIZE);
    Schur(fA, fV, SIZE);
    CArrayConvert(fV, V, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            Vt[i][j] = V[j][i];
        }
    }
    fVt = FArrayConvert(Vt, SIZE);
    
    fQ = FArrayConvert(Q, SIZE);
    fP = FArrayConvert(P, SIZE);
    fQP = FArrayConvert(QP, SIZE);
    MatMult(fQ, fP, fQP, SIZE);
    CArrayConvert(fQP, QP, SIZE);
    
    tr = Trace(QP, SIZE);
    //cout << "trace(QP) = " <<tr <<endl;
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            G[i][j] = QP[i][j]/tr;
        }
    }
    
    //  normalise by multiplying by V and V'(pre- and post-)
    fG = FArrayConvert(G, SIZE);
    fFoo = FArrayConvert(foo, SIZE);
    fGrad = FArrayConvert(gradMat, SIZE);
    MatMult(fV, fG, fFoo, SIZE);
    MatMult(fFoo, fVt, fGrad, SIZE);
    CArrayConvert(fGrad, gradMat, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] V[i];
        delete[] G[i];
        delete[] foo[i];
    }
    
    /***********************************************************************************/
    /***********************************************************************************/
    
    
    ofstream opfile3;
    opfile3.open("TESTQP.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile3 << QP[i][j] << "  ";
        }
        opfile3 <<endl;
    }
    
    ofstream opfile4;
    opfile4.open("TESTQ.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile4 << Q[i][j] << "  ";
        }
        opfile4 <<endl;
    }
    
    ofstream opfile5;
    opfile5.open("TESTP.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile5 << P[i][j] << "  ";
        }
        opfile5 <<endl;
    }
    
    //cout <<"exit of FormGradMat()" <<endl; 
    
    return;
}


/*  Performs gradient descent on W, based on gradients in gradMat  */
void GradDescent(double *W[], double *gradMat[], int SIZE)
{
    double descentRate;
    double wOld;

    // Determines how far down the gradient each iteration moves
    descentRate = 10;

    for(int i=0; i<SIZE; i++)
    {
        for(int j=(SIZE/2); j<SIZE; j++)
        {
            if(i != j)
            {
                wOld = W[i][j];
                W[i][j] -= descentRate*gradMat[i][j];
                if((i==80) && (j==76))
                {
                    //cout << "Wij updated from " <<wOld <<" to " <<W[i][j] <<endl;
                    //cout << "with gradient = " <<gradMat[i][j] <<endl;
                }
            }
        }
    }
    
    return;
}

