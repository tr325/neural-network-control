/*  Contains functions to generate an appropriate W matrix.  */

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<time.h>
#include"utilities.h"
#include"generateW.h"

using namespace std;

void GenerateWMat(double *W[], int *B[], int inhibCols, int SIZE)
{
    double inhibSparce;
    double exSparce;
    double inhibConst;
    double exConst; 
    int    exCols;
    double gamma;
    double exOmega;
    double inOmega;
    double specRad;
    double omSq;
    double exSq;
    double inSq;
    double sqrtN;
    double cMin;
    double kappa;

    srand(time(NULL));  //seeds the rand() function with the UNIX time

    while(inhibCols > SIZE)
    {
        cout << "Inhibitory columns must be less than SIZE, try again:" <<endl;
        cin >> inhibCols;
    }
    exCols = SIZE - inhibCols;
    inhibSparce = 0.2;
    exSparce = inhibSparce;  // allows the spectral radius formula used below (apparently easy to alter...)
    specRad = 0.4;
    kappa = 0.5;  // parameter for correllation antisymmetry

        
    // Only valid for constant sparcity.  Rajan for more detail (variance etc).
    sqrtN = double(sqrt(SIZE));
    omSq = specRad*specRad/(inhibSparce*(1-inhibSparce));
    exSq = omSq*inhibCols/exCols;
    inSq = omSq*exCols/inhibCols;
    exConst = sqrt(exSq)/sqrtN;
    inhibConst = -sqrt(inSq)/sqrtN;
    
    cMin = -exSparce/(1-exSparce);
    
    //~ gamma = 3.00;   // From Biology
    //~ 
    //~ exOmega = specRad/(sqrt(exSparce*(1-exSparce)*(1+gamma*gamma)/2));
    //~ inOmega = specRad/(sqrt(inhibSparce*(1-inhibSparce)*(1+gamma*gamma)/2));
    //~ exConst = exOmega/sqrt(SIZE);
    //~ inhibConst = -gamma*inOmega/sqrt(SIZE);
    //~ 
    
    //~ exConst = 1.054;  // What Guillaume sets it to in his paper
    //~ inhibConst = -((SIZE-inhibCols)*gamma*exSparce*exConst)/(inhibSparce*inhibCols);

    initB(B, cMin, exSparce, kappa, exCols, SIZE);

    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            if(j < exCols)
            {
                if(B[i][j])
                {
                    W[i][j] = exConst;
                }
                else
                {
                    W[i][j] = 0;
                }
            }
            else
            {
                if(B[i][j])
                {
                    W[i][j] = inhibConst;
                    B[i][j] = -1;
                }
                else
                {
                    W[i][j] = 0;
                }
            }
        }
    }

   /* 
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<exCols; j++)
        {
            if(randx() < exSparce)
            {
                W[i][j] = exConst;
            }
            else
            {
                W[i][j] = 0;
            }
        }
    }
    for(int i=0; i<SIZE; i++)
    {
        for(int j=exCols; j<SIZE; j++)
        {
            if(randx() < inhibSparce)
            {
                W[i][j] = inhibConst;
            }
            else
            {
                W[i][j] = 0;
            }
        }
    }
    */
    //EnforceDale(W, B, inhibCols, SIZE);

    OutputMat("generatedTestB.ascii", B, SIZE);
    
    return;
}

/*  Generate number between 0 and 1 */
double randx(void)
{
    double ret;

    ret = double(rand() % 1000)/1000;

    return ret;
}


/*  Initialises the B matrix, using pairwise correlation probabilities laid out in G's paper
 *  Only valid (in this form) for const sparcity (hence only one sparcity passed)
 *  skews towards positive symmetry                                               */
void initB(int *B[], double cMin, double exSparce, double kappa, int exCols, int SIZE)
{
    double recipTrue;
    double recipFalse;
    double cijSameNType;
    double cijOppNType;

    cijSameNType = kappa;      //cMax = 1
    cijOppNType = kappa*cMin;        

    //~ cout << "cij same neuron type = " << cijSameNType <<endl;
    //~ cout << "cij opp neuron type = " << cijOppNType <<endl;
    srand(time(NULL));
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            // Upper right triangle
            if(i < j)
            {
                if(randx() < exSparce)
                {
                    //~ cout << "upper right initilised " <<i << ", " <<j <<endl;
                    B[i][j] = 1;
                }
                else
                {
                    B[i][j] = 0;
                }
            }
            //  Lower left triangle
            else if(i > j)
            {
                // corresponding element set 
                if(B[j][i])
                {
                    // for E->E and I->I synapses
                    if(((i < exCols) && (j < exCols)) || ((i > exCols) && (j > exCols)))
                    {
                        recipTrue = exSparce + cijSameNType*(1 - exSparce);
                        if(randx() < recipTrue)
                        {
                            //~ cout << "lower left (same type) initilised " <<i << ", " <<j <<endl;
                            B[i][j] = 1;
                        }
                        else
                        {
                            B[i][j] = 0;
                        }
                    }
                    // for E->I and I->E synapses
                    else
                    {
                        recipTrue = exSparce + cijOppNType*(1 - exSparce);
                        if(randx() < recipTrue)
                        {
                            //~ cout << "lower left (opp type) initilised " <<i << ", " <<j <<endl;
                            B[i][j] = 1;
                        }
                        else
                        {
                            B[i][j] = 0;
                        }
                    }
                }
                // corresponding element not set
                else
                {
                                        // for E->E and I->I synapses
                    if(((i < exCols) && (j < exCols)) || ((i > exCols) && (j > exCols)))
                    {
                        recipFalse = exSparce*(1 - cijSameNType);
                        if(randx() < recipTrue)
                        {
                            B[i][j] = 1;
                        }
                        else
                        {
                            B[i][j] = 0;
                        }
                    }
                    // for E->I and I->E synapses
                    else
                    {
                        recipTrue = exSparce*(1 - cijOppNType);
                        if(randx() < recipTrue)
                        {
                            B[i][j] = 1;
                        }
                        else
                        {
                            B[i][j] = 0;
                        }
                    }
                }
            }
            // diagonal elements
            else
            {
                B[i][j] = 0;
            }
        }
    }

    return;
}

/*  Checks that the input matrix obeys Dale's Law and populates B matrix  
 *  NB: Should only be called once (at start).  The reparameterisation 
 *  should take care of keeping it enforced during the optimization.   */
void EnforceDale(double *W[], int *B[], int inhibNum, int SIZE)
{
    double exNum;
    exNum = SIZE - inhibNum;
    
    // Sort out the excitatory columns
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<exNum; j++)
        {
            if(W[i][j] <= 0)
            {
                W[i][j] = 0; 
                B[i][j] = 0;
            }
            else
            {
                B[i][j] = +1;
            }
        }
    }
    // Sort out inhibitory columns
    for(int i=0; i<SIZE; i++)
    {
        for(int j=exNum; j<SIZE; j++)
        {
            if(W[i][j] >= 0)
            {
                W[i][j] = 0;
                B[i][j] = 0;
            }
            else
            {
                B[i][j] = -1;
            }
        }
    }
    // Force diagonal elements to 0
    for(int i=0; i<SIZE; i++)
    {
        W[i][i] = 0;
        B[i][i] = 0;
    }    
        
    return;
}
