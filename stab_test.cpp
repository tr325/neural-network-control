//  Stabilization of W main file 

#include<iostream>
#include<iomanip>
#include<fstream>
#include"generateW.h"
#include"optimise.h"
#include"utilities.h"

using namespace std;


int main()
{
    int SIZE;
    double eps;
    int inhibCols;

    cout<< "Enter the dimension of matrix W:" <<endl;
    cin >> SIZE;
    cout << "How many inhibitory columns?" <<endl;
    cin >> inhibCols;
    eps = 0.01;

    cout << setprecision(7); 
    double *W[SIZE];
    int *B[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
        B[i] = new int[SIZE];
    }

    GenerateWMat(W, B, inhibCols, SIZE);
    OutputMat("generatedW.ascii", W, SIZE);
    OutputMat("generatedB.ascii", B, SIZE);
    OptimiseWMat(W, B, eps, inhibCols, SIZE);
    OutputMat("stabilizedW.ascii", W, SIZE);
    OutputMat("stabilizedB.ascii", B, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] W[i];
        delete[] B[i];
    }
    
    return 0;
}
