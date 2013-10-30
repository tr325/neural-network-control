// Test file fpor utilities.cpp

#include<iostream>
#include"utilities.h"

using namespace std; 

extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);

int main()
{
    double* p;
    double *A[3];
    double *B[3];


    int dimA = 3;
    int dimB = 3;
    
    for(int i=0; i<3; i++)
    {
        A[i] = new double[3]();
        B[i] = new double[3]();
    }
        
    //cout << A[1][1] << endl;
    
    p = lyap_P(A, dimA, B, dimB);
    
    
  /*  
    cout << p[0][0 <<endl;
    cout << p[1][0] <<endl;
    cout << p[2][0] <<endl; 
    */
    
    return 0; 
}
