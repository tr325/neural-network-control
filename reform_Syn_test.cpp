/*  Tests the ReformSyn() function          */

#include<iostream>
#include<iomanip>
#include<fstream>
#include"optimise.h"
#include"utilities.h"
#include"generateW.h"

using namespace std;


int main()
{
    ifstream ipfile;
    int SIZE;
    double eps;

    ipfile.open("refSynTestInput.ascii");
    SIZE = MatSize(ipfile);
    
    double *W[SIZE];
    double *V[SIZE];
    int *B[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
        V[i] = new double[SIZE];
        B[i] = new int[SIZE];
    }
    
    LoadMat(ipfile, W, SIZE);
    EnforceDale(W, B, SIZE, SIZE);
    Reparam(W, B, V, SIZE);
    ReformSyn(V, B, SIZE, SIZE);
    OutputMat("refSynTestV.ascii", V, SIZE);
    OutputMat("refSynTestB.ascii", B, SIZE);
    
    return 0;
}
    
    
