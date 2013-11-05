// Test file for utilities.cpp

#include<fstream>
#include"utilities.h"
#include<iostream>

using namespace std; 

extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);


int main()
{
  
    double *test[SIZE], *test2[SIZE], *test3[SIZE];
    double epsilon = 1;
    
    
    for(int i=0; i<SIZE; i++)
    {
        test[i] = new double[SIZE];
        test2[i] = new double[SIZE];
        test3[i] = new double[SIZE]; 
        for(int j=0; j<SIZE; j++)
        {
            test[i][i] = 1;
            test2[i][i] = 1;
            test3[i][j] = 1; 
        }
    }
    test[1][2] = 4; 
    test[2][1] = 3;
    test2[2][1] = 10;
    
    cout << "Epsilon = " << epsilon << endl; 

    
    ifstream myfile;
    myfile.open("inputMatrix.ascii");
    loadMat(myfile, test, SIZE);
    double ssa = SimpleSSA(test, epsilon);    
    
    return 0; 
}
