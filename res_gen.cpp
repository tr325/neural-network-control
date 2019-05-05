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
    double *Q[SIZE];
    double *P[SIZE];
    double *Wref[SIZE];
    double *A[SIZE];
    double *sV[SIZE];
    int *I[SIZE];
    int *B[SIZE];
    int *Bref[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        W[i] = new double[SIZE];
        Wref[i] = new double[SIZE];
        B[i] = new int[SIZE];
        Bref[i] = new int[SIZE];
        A[i] = new double[SIZE];
        I[i] = new int[SIZE];
        sV[i] = new double[SIZE];
        P[i] = new double[SIZE];
        Q[i] = new double[SIZE];
    }

    /*************** Data for smoothness of SA plot **********/
    double dimVar;
    double sA;
    double ssA1;
    double ssA2;
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            if(i == j)
            {
                I[i][j] = 1;
            }
            else
            {
                I[i][j] = 0;
            }
        }
    }

    GenerateWMat(W, B, inhibCols, SIZE, false, true, 0);

    ofstream ssafile;
    ssafile.open("smoothness.ascii");
    
    for(dimVar = -50; dimVar < 50; dimVar = dimVar + 0.1)
    {
        W[18][5] = dimVar;
        
        ssA1 = SimpleSSA(W, P, Q, sV, 0.001, SIZE);
        ssA2 = SimpleSSA(W, P, Q, sV, 0.01, SIZE);

        for(int i=0; i<SIZE; i++)
        {
            for(int j=0; j<SIZE; j++)
            {
                A[i][j] = W[i][j] - I[i][j];
            }
        }

        Schur(A, sV, SIZE);
        sA = MaxDiag(A, SIZE);


        cout << "sa = " << sA << " ssA = " <<ssA1 <<endl;
        ssafile << sA << "   " << ssA1 << "   " <<ssA2 <<endl;
    }

    ssafile.close();
//~ 
    //~ 
    
    /************************************************
     * change the numbers (on input), but CHANGE FILENAME TOO
     * ***********************************************/
    
    //~ stringstream ss;
    //~ string fname1 = "limW";  //edit if you do a different run to avoid data loss
    //~ string fname3 = "newBeta";
    //~ string fname2 = ".ascii";
    //~ string ref = "REF";
    //~ string fname;
    //~ string reffname;
    //~ 
    //~ double beta;
    //~ bool wholeNet = false;
    //~ bool pairwise = false;
    //~ bool wPlast = true;
    //~ char *opfile;
    //~ double kappa;
    //~ double betas[5] = {1, 1.5, 2, 2.5, 3};
    //~ 
//~ 
    //~ 
    //~ GenerateWMat(W, B, inhibCols, SIZE, pairwise, true, kappa);
    //~ ss << fname1 <<ref << fname2 <<fname3;
    //~ reffname = ss.str();
    //~ opfile = &reffname[0];
    //~ OutputMat(opfile, W, SIZE);
    //~ ss.str("");

    //~ for(int i=0; i<SIZE; i++)
    //~ {
        //~ for(int j=0; j<SIZE; j++)
        //~ {
            //~ Wref[i][j] = W[i][j];
            //~ Bref[i][j] = B[i][j];
        //~ }
    //~ }
    //~ 
    //~ for(int k=0; k<5; k++)
    //~ {
//~ 
        //~ kappa = 0.2*k;    
        //~ GenerateWMat(W, B, inhibCols, SIZE, pairwise, wholeNet, kappa);
        //~ beta = betas[k];
        //~ cout << "******************************************" <<endl;
        //~ cout << "BETA= " <<beta << endl;
        //~ 
        
        //~ 
        //~ for(int n=0; n<25; n++)
        //~ {
            //~ // start from ref network again for fixed E fixed I trials
            //~ // such as beta tests (scaling network)
            //~ for(int i=0; i<SIZE; i++)
            //~ {
                //~ for(int j=0; j<SIZE; j++)
                //~ {
                    //~ W[i][j] = beta*Wref[i][j];
                    //~ B[i][j] = Bref[i][j];
                //~ }
            //~ }
            //~ //~
            //~ 
            //~ cout << endl <<endl << "OPTIMISATION NUMBER :" <<n+1 <<endl;
            //~ ss.str("");
            //~ GenerateWMat(W, B, inhibCols, SIZE, pairwise, wholeNet, kappa);
            //~ OptimiseWMat(W, B, eps, inhibCols, wPlast, SIZE);
            //~ ss << fname1;
            //~ ss << n;
            //~ ss << fname3;
            //~ ss << beta;
            //~ ss << fname2;
            //~ fname = ss.str();
            //~ opfile = &fname[0];
            //~ OutputMat(opfile, W, SIZE);
        //~ }
    //~ }

    for(int i=0; i<SIZE; i++)
    {
        delete[] W[i];
        delete[] B[i];
    }
    
    return 0;
}
