// Test file for reparameterisation and W recalculation

#include<iostream>
#include<iomanip>
#include<fstream>
#include"optimise.h"
#include"utilities.h"

using namespace std;

int main()
{
    int SIZE;
    ifstream ipfile; 
    ipfile.open("repTestInput.ascii");
    SIZE = MatSize(ipfile);
    
    double  *W[SIZE];
    double  *Wnew[SIZE];
    double  *V[SIZE];
    int     *B[SIZE];
    int     inhibNum;
    
    for(int i =0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
        Wnew[i] = new double[SIZE];
        V[i] = new double[SIZE];
        B[i] = new int[SIZE];
    }
    
    LoadMat(ipfile, W, SIZE);
    
    inhibNum = 10;
    EnforceDale(W, B, inhibNum, SIZE);
    Reparam(W, B, V, SIZE);
    RecalcW(Wnew, B, V, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            if( (Wnew[i][j] - W[i][j])  >  0.0000001)
            {
                double e = (Wnew[i][j] - W[i][j]);
                cout << "Error - out by " <<e <<endl;
            }
        }
    }
    
    OutputMat("reptestW.ascii", W, SIZE);
    OutputMat("reptestWnew.ascii", Wnew, SIZE);
    
    
    return 0;
}
    
