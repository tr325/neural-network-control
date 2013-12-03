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
    double *fW;
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
    /*
    fW = FArrayConvert(W, SIZE);
    Schur(fW, SIZE);
    CArrayConvert(fW, W, SIZE);
    */
    ssaFIXED = SimpleSSA(W, P, Q, eps, SIZE); 
    cout << "Original ssa value = " <<ssaFIXED <<endl;
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            A[i][j] = W[i][j] - ssaFIXED*I[i][j];
        }
    }
    
    ofstream opfile4; 
    opfile4.open("TESTA.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile4 << A[i][j] <<"  ";
        }
        opfile4 <<endl;
    }
    opfile4.close();
    
    ofstream opfile3;
    opfile3.open("TESTPssa.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile3 << P[i][j] <<"  ";
        }
        opfile3 <<endl;
    }
    opfile3.close();
    
    cout << "*******************************************" <<endl;
    FormGradMat(gradMat, P, Q, A, QP, ssaFIXED, SIZE);

    ofstream opfile2;
    opfile2.open("TESTGradMat.ascii");
    opfile2.precision(10);
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile2 << gradMat[i][j] << "  ";
        }
        opfile2 << endl;
    }
    opfile2.close();
    
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
    opfile.precision(10);    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile << numGradMat[i][j] << "  ";
        }
        opfile << endl; 
    }
    opfile.close();
       
    cout << "true 1 = " << gradMat[0][0] <<endl;
    cout << "numerical 1 = " << numGradMat[0][0] <<endl;
    cout << "true 2 = " << gradMat[2][2] <<endl;
    cout << "numerical 2 = " << numGradMat[2][2] <<endl;

    
    return 0;
}
