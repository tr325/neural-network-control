/*  Calculates the gradient matrix by numerical methods */

#include<fstream>
#include<iostream>
#include"utilities.h"
#include"optimise.h"


using namespace std;

int main()
{
    ifstream ipfile;
    int SIZE;
    

    ipfile.open("repTestInput.ascii");
    SIZE = MatSize(ipfile);
    
    cout << "SIZE = " <<SIZE <<endl;
    
    double *W[SIZE];
    double *P[SIZE];
    double *Q[SIZE];
    double *A[SIZE];
    double *I[SIZE];
    double *gradMat[SIZE];
    double *numGradMat[SIZE];
    double *QP[SIZE];
    int *B[SIZE];
    double eps;
    double delta;
    double ssa;
    double ssaFIXED;
    
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
        A[i] = new double[SIZE];
        I[i] = new double[SIZE];
        P[i] = new double[SIZE];
        Q[i] = new double[SIZE];
        QP[i] = new double[SIZE];
        B[i] = new int[SIZE];
        gradMat[i] = new double[SIZE];
        numGradMat[i] = new double[SIZE];
    }
    
    for(int i=0; i<SIZE; i++)
    {
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
    
    LoadMat(ipfile, W, SIZE);
    delta = 0.00001; 
    eps = 0.001;
    
    EnforceDale(W, B, 10, SIZE);
    ssaFIXED = SimpleSSA(W, P, Q, eps, SIZE); 
    cout << "Original ssa value = " <<ssaFIXED <<endl;
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            A[i][j] = W[i][j] - ssaFIXED*I[i][j];
        }
    }
    
    MatMult(Q, P, QP, SIZE);
    cout << "*******************************************" <<endl;
    FormGradMat(gradMat, W, A, QP, ssaFIXED, SIZE);

    
    //cout << "Forming numerical gradMat..." << endl;
    //  Perturbs each element of W by a small amount delta, recalculates ssa
    //  Gradient is then (change in ssa)/delta
    
    /*
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            W[i][j] += delta;
            ssa = SimpleSSA(W, P, Q, eps, SIZE);
            numGradMat[i][j] = (ssa-ssaFIXED)/(delta);
            W[i][j] -= delta;           
        }
        cout << "Completed line " << i <<endl; 
    }    
    */
    /*   
    cout << "true 1 = " << gradMat[0][0] <<endl;
    cout << "numerical 1 = " << numGradMat[0][0] <<endl;
    cout << "true 2 = " << gradMat[2][2] <<endl;
    cout << "numerical 2 = " << numGradMat[2][2] <<endl;
    */
    
    return 0;
}
