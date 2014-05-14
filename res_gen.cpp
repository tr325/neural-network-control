// Generation of results main file.  Commented sections give different results

#include<iostream>
#include<iomanip>
#include<fstream>
#include<stdlib.h>
#include<string>
#include<sstream>
#include"generateW.h"
#include"optimise.h"
#include"utilities.h"

using namespace std;


int main(int argc, char* argv[])
{
    double eps;
    
    // decode arguments
    if(argc != 3) {
        printf("You must provide 2 arguments - dimension of W and the number of inhibitory columns\n");
        exit(0);
    }
    // report settings
    printf("SIZE = %i\n", atoi(argv[1]));
    printf("inhibCols = %i\n", atoi(argv[2]));
    
    int SIZE = atoi(argv[1]);
    int inhibCols = atoi(argv[2]);
    
    eps = 0.001;

    cout << setprecision(7); 
    double *W[SIZE];
    int *B[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
        B[i] = new int[SIZE];
    }

    /************************************************
     * Fixed E network, different I networks, no plasticity
     * (change the numbers if you want (on input), but CHANGE FILENAME TOO
     * ***********************************************/
    
    stringstream ss;
    string fname1 = "fixedE8020noplast";
    string fname2 = ".ascii";
    string fname;
    bool wholeNet = false;
    bool pairwise = false;
    bool wPlast = false;
    char *opfile;
    GenerateWMat(W, B, inhibCols, SIZE, pairwise, true);
    for(int i=0; i<25; i++)
    {
        ss.str("");
        GenerateWMat(W, B, inhibCols, SIZE, pairwise, wholeNet);
        OptimiseWMat(W, B, eps, inhibCols, wPlast, SIZE);
        ss << fname1;
        ss << i;
        ss << fname2;
        fname = ss.str();
        opfile = &fname[0];
        OutputMat(opfile, W, SIZE);
    }
    

    /**************************************************
     * Fixed E network, different initial I nets, with plasticity
     * *********************************************************/
    /*
    stringstream ss;
    string fname1 = "fixedE5050withPlast";
    string fname2 = ".ascii";
    string fname;
    bool wholeNet = false;
    bool pairwise = false;
    bool wPlast = true;
    char *opfile;
    GenerateWMat(W, B, inhibCols, SIZE, pairwise, true);
    for(int i=0; i<25; i++)
    {
        ss.str("");
        GenerateWMat(W, B, inhibCols, SIZE, pairwise, wholeNet);
        OptimiseWMat(W, B, eps, inhibCols, wPlast, SIZE);
        ss << fname1;
        ss << i;
        ss << fname2;
        fname = ss.str();
        opfile = &fname[0];
        OutputMat(opfile, W, SIZE);
    }
    */

    for(int i=0; i<SIZE; i++)
    {
        delete[] W[i];
        delete[] B[i];
    }
    
    return 0;
}
