// Test file for utilities.cpp

#include<fstream>
#include"utilities.h"
#include<iostream>

using namespace std; 

extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);


int main()
{
  
    double *test[SIZE], *test2[SIZE];
    double epsilon = 0.01;
    ofstream opfile;
    
    for(int i=0; i<SIZE; i++)
    {
        test[i] = new double[SIZE];
        test2[i] = new double[SIZE];
    }
   

    cout << "Epsilon = " << epsilon << endl;
    cout << "SIZE = " <<SIZE <<endl;

    ifstream myfile;
    myfile.open("inputMatrix.ascii");
    //myfile.open("inputIdent.ascii");
    //myfile.open("TESTinputMatrix.ascii");
    loadMat(myfile, test, SIZE);
    /*
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            test2[i][j] = test[j][i];
        }
    }
    
    double *fTest, *fTest2, *fLyapResult; 
    fTest2 = FArrayConvert(test2, SIZE);
    fTest = FArrayConvert(test, SIZE);
    Schur(fTest);
    Schur(fTest2);
    fLyapResult = Lyap(fTest, fTest2);
    CArrayConvert(fLyapResult, test2, SIZE);
    CArrayConvert(fTest, test, SIZE);
    
    opfile.open("SchurInputMatrix");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile << test[i][j] << "  ";
        }
        opfile <<endl;
    }
    
    
    opfile.open("LyapTestResult");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile << test2[i][j] << "  ";
        }
        opfile <<endl;
    }
    */
    
    //opfile.open("epsVsSSA.ascii", ios::app);
    double ssa;
    //for(epsilon = 0.2; epsilon < 1; epsilon += 0.1)
    {
        cout << "Epsilon = " <<epsilon << endl; 
        ssa = SimpleSSA(test, epsilon);    
        //opfile << ssa << " " << epsilon << endl; 
    }    
    
    return 0; 
}
