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
    

    ipfile.open("TESTinputMatrix.ascii");
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
    double *fP;
    double *fQ;
    double *fQP;
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
    
    loadMat(ipfile, W, SIZE);
    delta = 0.00001; 
    eps = 0.001;
    
    ssaFIXED = SimpleSSA(W, P, Q, eps, SIZE); 
    cout << "Original ssa value = " <<ssaFIXED <<endl;
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            A[i][j] = W[i][j] - ssaFIXED*I[i][j];
        }
    }

    fP = FArrayConvert(P, SIZE);
    fQ = FArrayConvert(Q, SIZE);
    fQP = FArrayConvert(QP, SIZE);
    MatMult(fQ, fP, fQP, SIZE);
    delete[] fP;
    delete[] fQ;
    CArrayConvert(fQP, QP, SIZE);
    
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
    opfile4.open("TESTP.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile4 << P[i][j] << "  ";
        }
        opfile4 <<endl;
    }
    ofstream opfile5;
    opfile5.open("TESTQ.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile5 << Q[i][j] << "  ";
        }
        opfile5 <<endl;
    }
    
    delete[] fQP;
    FormGradMat(gradMat, P, Q, A, QP, ssaFIXED, SIZE);

    ofstream opfile2;
    opfile2.open("TESTGradMat.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile2 << gradMat[i][j] << "  ";
        }
        opfile2 << endl;
    }
    
    cout << "Forming numerical gradMat..." << endl;
    //  Perturbs each element of W by a small amount delta, recalculates ssa
    //  Gradient is then (change in ssa)/delta
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
    
    ofstream opfile; 
    opfile.open("TESTNumGradMat.ascii");
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile << numGradMat[i][j] << "  ";
        }
        opfile << endl; 
    }
       
    cout << "true 1 = " << gradMat[15][7] <<endl;
    cout << "estimated 1 = " << numGradMat[15][7] <<endl;
    cout << "ratio 1 = " << gradMat[2][2]/numGradMat[2][2] <<endl;
    cout << "ratio 2 = " <<gradMat[15][7]/numGradMat[15][7] <<endl;
    cout << "difference 1 = " << gradMat[2][2]-numGradMat[2][2] <<endl;
    cout << "difference 2 = " <<gradMat[15][7]-numGradMat[15][7] <<endl;
    
    return 0;
}
