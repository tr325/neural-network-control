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
    double *Wref[SIZE];
    int *B[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
        Wref[i] = new double[SIZE];
        B[i] = new int[SIZE];
    }

    /************************************************
     * change the numbers (on input), but CHANGE FILENAME TOO
     * ***********************************************/
    
    stringstream ss;
    string fname1 = "fEfI8020withplast";  //edit if you do a different run to avoid data loss
    string fname3 = "Kappa";
    string fname2 = ".ascii";
    string ref = "REF";
    string fname;
    string reffname;
    int beta;
    bool wholeNet = true;
    bool pairwise = true;
    bool wPlast = true;
    char *opfile;
    double kappa;

    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            Wref[i][j] = W[i][j];
        }
    }

    ss << fname1 <<ref << fname2 <<fname3;
    reffname = ss.str();
    opfile = &reffname[0];
    OutputMat(opfile, W, SIZE);
    ss.str("");

    for(int k=0; k<6; k++)
    {

        kappa = 0.2*k;    
        GenerateWMat(W, B, inhibCols, SIZE, pairwise, wholeNet, kappa);
        
        for(int n=0; n<25; n++)
        {
            //~ // start from ref network again for fixed E fixed I trials
            //~ for(int i=0; i<SIZE; i++)
            //~ {
                //~ for(int j=0; j<SIZE; j++)
                //~ {
                    //~ W[i][j] = beta*Wref[i][j];
                //~ }
            //~ }
            //~ 
            cout << endl <<endl << "OPTIMISATION NUMBER :" <<n+1 <<endl;
            ss.str("");
            GenerateWMat(W, B, inhibCols, SIZE, pairwise, wholeNet, kappa);
            OptimiseWMat(W, B, eps, inhibCols, wPlast, SIZE);
            ss << fname1;
            ss << n;
            ss << fname3;
            ss << kappa;
            ss << fname2;
            fname = ss.str();
            opfile = &fname[0];
            OutputMat(opfile, W, SIZE);
        }
    }

    for(int i=0; i<SIZE; i++)
    {
        delete[] W[i];
        delete[] B[i];
    }
    
    return 0;
}
