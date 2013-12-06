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
    double *fA;
    bool conv;
    double sa;
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
    descentRate = 10;   // determines how far down slope each iteration updates
    
    ofstream resFile;
    resFile.open("SsaOptimise.ascii");
    
    // loops until the gradient descent algorithm has converged
    while(!conv)
    {
        //cout << "Optimisation loopcount = " <<loopcount << endl; 
        // solve SSA for current W matrix
        ssaOld = ssa;
        ssa = SimpleSSA(W, P, Q, eps, SIZE);

        for(int i=0; i<SIZE; i++)
        {
            for(int j=0; j<SIZE; j++)
            {
                A[i][j] = W[i][j] - ssa*I[i][j];
            }
        }
        
        // find the gradient matrix, d(SSA)/dW;
        MatMult(Q, P, QP, SIZE);
        FormGradMat(gradMat, A, QP, ssa, SIZE);
        
        // perform gradient descent to optimise W
        GradDescent(W, gradMat, descentRate, SIZE);
        
        if(abs(ssa - ssaOld) < precision)
        {
            cout << "Converged" <<endl;
            conv = true;
        }
        /*
        if(loopcount == 1000)
        {
            cout <<"That's enough i think..." <<endl;
            conv = true;
        }
        */
        
        resFile << ssa <<endl; 
        
        cout << "SSA on loop " << loopcount << ": " <<ssa <<endl; 
        loopcount++;
    }
    resFile.close();

    for(int i=0; i<SIZE; i++)
    {
        delete[] P[i];
        delete[] Q[i];
        delete[] gradMat[i];
        delete[] QP[i];
        delete[] A[i];
        delete[] I[i];
    }
    delete[] fP;
    delete[] fQ;
    delete[] fQP;
    delete[] fA;

    return;
}


/*  Forms the gradient matrix           */
void FormGradMat(double *gradMat[], double *A[], double *QP[], double ssa, int SIZE)
{
    double *V[SIZE];
    double *Vt[SIZE];
    double *G[SIZE];
    double *foo[SIZE];
    double tr;
    
    for(int i=0; i<SIZE; i++)
    {
        V[i] = new double[SIZE];
        Vt[i] = new double[SIZE];
        G[i] = new double[SIZE];
        foo[i] = new double[SIZE];
    }
    
    //  Calculate the matrix of schur vectors, V
    Schur(A, V, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            Vt[i][j] = V[j][i];
        }
    }
    
    tr = Trace(QP, SIZE);    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            G[i][j] = QP[i][j]/tr;
        }
    }
    
    //  normalise by multiplying by V and V'(pre- and post-)
    MatMult(V, G, foo, SIZE);
    MatMult(foo, Vt, gradMat, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] V[i];
        delete[] G[i];
        delete[] foo[i];
    }
    

    /***********************************************************************************/
    /***********************************************************************************/
    
    /*
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
    opfile3.close();
    
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
    opfile4.close();
    
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
    opfile5.close()
    */
    //cout <<"exit of FormGradMat()" <<endl; 
    
    return;
}


/*  Performs gradient descent on W, based on gradients in gradMat  */
void GradDescent(double *W[], double *gradMat[], double descentRate, int SIZE)
{
    for(int i=0; i<SIZE; i++)
    {
        for(int j=(SIZE/2); j<SIZE; j++)
        {
            if(i != j)
            {
                W[i][j] -= descentRate*gradMat[i][j];
            }
        }
    }
    
    return;
}

