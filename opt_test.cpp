//  Test file for optimisation functions

#include<iostream>
#include<iomanip>
#include<fstream>
#include"optimise.h"
#include"utilities.h"

using namespace std;


int main()
{
    ifstream ipfile;
    int SIZE;
    double eps;

    ipfile.open("gradOpt.ascii");
    SIZE = MatSize(ipfile);
    eps = 0.01;

    cout << setprecision(7); 
    double *W[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
    }

    LoadMat(ipfile, W, SIZE);

    OptimiseWMat(W, eps, SIZE);     
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] W[i];
    }
    
    return 0;
}
