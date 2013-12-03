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
    double *fW;

    ipfile.open("gradOpt.ascii");
    SIZE = MatSize(ipfile);
    eps = 0.01;

    cout << setprecision(7); 
    double *W[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
    }

    loadMat(ipfile, W, SIZE);
    /*
    fW = FArrayConvert(W, SIZE);
    Schur(fW, SIZE);
    CArrayConvert(fW, W, SIZE);
    */
    OptimiseWMat(W, eps, SIZE);     


    ofstream opfile;
    opfile.precision(7);
    opfile.open("editedW.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile << W[i][j] << "  ";
        }
        opfile <<endl;
    }
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] W[i];
    }
    
    return 0;
}
