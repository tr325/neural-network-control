//  Stabilization of W main file 

#include<iostream>
#include<iomanip>
#include<fstream>
#include<stdlib.h>
#include"generateW.h"
#include"optimise.h"
#include"utilities.h"

using namespace std;


int main(int argc, char* argv[])
{
    double eps;
    
    // decode arguments
    if(argc < 2) {
        printf("You must provide 2 arguments - dimension of W and the number of inhibitory columns\n");
        exit(0);
    }
    // report settings
    printf("SIZE = %i\n", atoi(argv[1]));
    printf("inhibCols = %i\n", atoi(argv[2]));
    
    int SIZE = atoi(argv[1]);
    int inhibCols = atoi(argv[2]);
    
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
