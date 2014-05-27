/*  Tests the GenerateW() function          */

#include<iostream>
#include<iomanip>
#include<fstream>
#include<stdlib.h>
#include"generateW.h"
#include"utilities.h"


using namespace std;


int main(int argc, char* argv[])
{
    // decode arguments
    if(argc != 3) {
        printf("You must provide two arguments - dimension of W and the number of inhibitory columns\n");
        printf("argc = %i\n", argc);
        return(0);
    }
    // report settings
    //~ printf("SIZE = %i\n", atoi(argv[1]));
    //~ printf("inhibCols = %i\n", atoi(argv[2]));
    
    int SIZE = atoi(argv[1]);
    int inhibCols = atoi(argv[2]);
    bool pairwise = false;
    bool wholeNet = true;
    double kappa = 1;
    
    double *W[SIZE];
    int     *B[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
        B[i] = new int[SIZE];
    }
    
    GenerateWMat(W, B, inhibCols, SIZE, pairwise, wholeNet, kappa);
    OutputMat("generatedTestW.ascii", W, SIZE);
    
    
    return 0;
}
