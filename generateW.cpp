/*  Contains functions to generate an appropriate W matrix.  */

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<time.h>
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
    double randx;
    double exOmega;
    double inOmega;
    double specRad;
    double omSq;
    double exSq;
    double inSq;
    double sqrtN;
    
    while(inhibCols > SIZE)
    {
        cout << "Inhibitory columns must be less than SIZE, try again:" <<endl;
        cin >> inhibCols;
    }
    exCols = SIZE - inhibCols;
    inhibSparce = 0.2;
    exSparce = inhibSparce;  // allows the spectral radius formula used below (apparently easy to alter...)
    specRad = 3.00;

    // Only valid for constant sparcity.  Rajan for more detail (variance etc).
    sqrtN = double(sqrt(SIZE));
    omSq = specRad*specRad/(inhibSparce*(1-inhibSparce));
    exSq = omSq*inhibCols/exCols;
    inSq = omSq*exCols/inhibCols;
    exConst = sqrt(exSq)/sqrtN;
    inhibConst = -sqrt(inSq)/sqrtN;
    
    
    
    //~ gamma = 3.00;   // From Biology
    //~ 
    //~ exOmega = specRad/(sqrt(exSparce*(1-exSparce)*(1+gamma*gamma)/2));
    //~ inOmega = specRad/(sqrt(inhibSparce*(1-inhibSparce)*(1+gamma*gamma)/2));
    //~ exConst = exOmega/sqrt(SIZE);
    //~ inhibConst = -gamma*inOmega/sqrt(SIZE);
    //~ 
    
    //~ exConst = 1.054;  // What Guillaume sets it to in his paper
    //~ inhibConst = -((SIZE-inhibCols)*gamma*exSparce*exConst)/(inhibSparce*inhibCols);
    srand(time(NULL));  //seeds the rand() function with the UNIX time
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<exCols; j++)
        {
            //Generates random number between 0 and 1
            randx = double(rand() % 1000)/1000;
            if(randx < exSparce)
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
            //Generates random number between 0 and 1
            randx = double(rand() % 1000)/1000;  
            if(randx < inhibSparce)
            {
                W[i][j] = inhibConst;
            }
            else
            {
                W[i][j] = 0;
            }
        }
    }
    
    EnforceDale(W, B, inhibCols, SIZE);
    
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
