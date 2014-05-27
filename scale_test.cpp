//  checking scaling of the W matrix stabilisation

#include<iostream>
#include<iomanip>
#include<fstream>
#include<stdlib.h>
#include"generateW.h"
#include"optimise.h"
#include"utilities.h"

using namespace std;


int main(void)
{
    double eps;
    
    int SIZE = 100;
    int inhibCols = 20;
    
    eps = 0.001;

    cout << setprecision(7); 
    double *W1[SIZE];
    double *W2[SIZE];
    double *W4[SIZE];
    double *W8[SIZE];
    double *W16[SIZE];
    double *Wref[SIZE];
    int *B[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        W1[i] = new double[SIZE];
        W2[i] = new double[SIZE];
        W4[i] = new double[SIZE];
        W8[i] = new double[SIZE];
        W16[i] = new double[SIZE];
        B[i] = new int[SIZE];
        Wref[i] = new double[SIZE];
    }

    char* s = "Wref.ascii";
    char t;
    char* u;
    char* stu;
    GenerateWMat(Wref, B, inhibCols, SIZE);
    OutputMat(s, Wref, SIZE);

    for(int j=0; j<SIZE; j++)
    {
        for(int k=0; k<SIZE; k++)
        {
            W1[j][k] = Wref[j][k];
            W2[j][k] = 2*Wref[j][k];
            W4[j][k] = 4*Wref[j][k];
            W8[j][k] = 8*Wref[j][k];
            W16[j][k] = 16*Wref[j][k];
        }
    }

    OptimiseWMat(W1, B, eps, inhibCols, SIZE);        
    OptimiseWMat(W2, B, eps, inhibCols, SIZE);        
    OptimiseWMat(W4, B, eps, inhibCols, SIZE);        
    OptimiseWMat(W8, B, eps, inhibCols, SIZE);        
    
    OutputMat("Wbeta1.ascii", W1, SIZE);
    OutputMat("Wbeta2.ascii", W2, SIZE);
    OutputMat("Wbeta4.ascii", W4, SIZE);
    OutputMat("Wbeta8.ascii", W8, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] W1[i];
        delete[] W2[i];
        delete[] W4[i];
        delete[] W8[i];
        delete[] W16[i];
        delete[] Wref[i];
        delete[] B[i];
    }
    
    return 0;
}
