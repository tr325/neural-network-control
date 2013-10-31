// Test file for utilities.cpp

#include<iostream>
#include"utilities.h"

using namespace std; 

extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);


int main()
{
  
    double *test[SIZE];
    double epsilon = 0.50;
    
    
    for(int i=0; i<SIZE; i++)
    {
        test[i] = new double[SIZE];
        for(int j=0; j<SIZE; j++)
        {
            test[i][i] = 1;
        }
    }
    test[1][2] = -4; 
    test[2][1] = 3;
    
    cout << "Epsilon = " << epsilon << endl; 
    cout << "test matrix:" << endl; 
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            cout << test[i][j] << "    ";
        }
        cout << endl; 
    }
    cout << endl;
    
    double ssa = SimpleSSA(test, epsilon);  
    
    return 0; 
}
