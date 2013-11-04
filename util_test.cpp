// Test file for utilities.cpp

#include<iostream>
#include"utilities.h"

using namespace std; 

extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);


int main()
{
  
    double *test[SIZE], *test2[SIZE], *test3[SIZE];
    double epsilon = 200;
    
    
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
    
    double *fTest, *fTest2, *fTest3;
    fTest = FArrayConvert(test, SIZE);
    fTest2 = FArrayConvert(test2, SIZE);
    fTest3 = FArrayConvert(test3, SIZE); 
    MatMult(fTest, fTest2, fTest3);
    CArrayConvert(fTest3, test3, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            cout << test3[i][j] << "    ";
        }
        cout << endl; 
    }
    cout << endl;
    double ssa = SimpleSSA(test, epsilon);  
    
    return 0; 
}
