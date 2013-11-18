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
    
    double *W[SIZE];
    double *P[SIZE];
    double *Q[SIZE];
    double *gradMat[SIZE];
    double *numGradMat[SIZE];
    double *PQ[SIZE];
    double *fP;
    double *fQ;
    double *fPQ;
    double eps;
    double delta;
    double ssa;
    
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
        P[i] = new double[SIZE];
        Q[i] = new double[SIZE];
        PQ[i] = new double[SIZE];
        gradMat[i] = new double[SIZE];
        numGradMat[i] = new double[SIZE];
    }

    loadMat(ipfile, W, SIZE);
    delta = 0.00001; 
    eps = 0.01;
    
        

    ssa = SimpleSSA(W, P, Q, eps, SIZE); 

    fP = FArrayConvert(P, SIZE);
    fQ = FArrayConvert(Q, SIZE);
    fPQ = FArrayConvert(PQ, SIZE);
    MatMult(fP, fQ, fPQ, SIZE);
    
    cout << "here" << endl;

    delete[] fP;
    delete[] fQ;
    CArrayConvert(fPQ, PQ, SIZE);
    delete[] fPQ;
    FormGradMat(gradMat, PQ, SIZE);
    
    ofstream opfile2;
    opfile2.open("TESTGradMat");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile2 << gradMat[i][j] << "  ";
        }
        opfile2 << endl;
    }
    
    cout << "Forming numerical gradMat..." << endl;
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            W[i][j] += delta;
            ssa = SimpleSSA(W, P, Q, eps, SIZE); 
            W[i][j] -= 2*delta; 
            ssa -= SimpleSSA(W, P, Q, eps, SIZE);
            numGradMat[i][j] = ssa/(2*delta);
            W[i][j] += delta;           
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
    
    return 0;
}
