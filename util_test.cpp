// Test file for utilities.cpp

#include<fstream>
#include"utilities.h"
#include<iostream>

using namespace std;
extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);


int main()
{
    ifstream ipfile;
    int SIZE;
    
    ipfile.open("inputMatrix.ascii");
    SIZE = MatSize(ipfile); 

    double epsilon = 0.02;
    ofstream opfile;
    
    double *test[SIZE];
    double *P[SIZE];
    double *Q[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        test[i] = new double[SIZE];
        P[i] = new double[SIZE];
        Q[i] = new double[SIZE];
    }
   
    cout << "Epsilon = " << epsilon << endl;
    cout << "SIZE = " <<SIZE <<endl;

    loadMat(ipfile, test, SIZE);
    
    opfile.open("epsVsSSA.ascii", ios::trunc);
    double ssa;
    for(epsilon = 0.001; epsilon < 0.1; epsilon += 0.001)
    {
        cout << "Epsilon = " <<epsilon << endl; 
        ssa = SimpleSSA(test, P, Q, epsilon, SIZE);    
        opfile << epsilon << "  " << ssa << endl; 
    }
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] test[i];
    }
    
    return 0; 
}
