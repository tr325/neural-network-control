// Test file for utilities.cpp

#include<fstream>
#include"utilities.h"
#include<iostream>

using namespace std;
extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);


int main()
{
    double *test[SIZE];
    double epsilon = 0.02;
    ofstream opfile;
    ifstream ipfile;
    
    for(int i=0; i<SIZE; i++)
    {
        test[i] = new double[SIZE];
    }
   

    cout << "Epsilon = " << epsilon << endl;
    cout << "SIZE = " <<SIZE <<endl;


    ipfile.open("inputMatrix.ascii");
    //myfile.open("inputIdent.ascii");
    //myfile.open("TESTinputMatrix.ascii");
    loadMat(ipfile, test, SIZE);
    
    opfile.open("epsVsSSA.ascii", ios::trunc);
    double ssa;
    for(epsilon = 0.001; epsilon < 0.1; epsilon += 0.001)
    {
        cout << "Epsilon = " <<epsilon << endl; 
        ssa = SimpleSSA(test, epsilon);    
        opfile << epsilon << "  " << ssa << endl; 
    }
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] test[i];
    }
    
    return 0; 
}
