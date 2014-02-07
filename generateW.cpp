/*  Contains functions to generate an appropriate W matrix.  */

#include<iostream>
#include<cmath>
#include<cstdlib>
#include"generateW.h"

using namespace std;

void GenerateWMat(double *W[], int *B[], int inhibCols, int SIZE)
{
    double inhibSparce;
    double exSparce;
    double inhibConst;
    double exConst; 
    int exCols;
    double gamma;
    double randx;
    
    while(inhibCols > SIZE)
    {
        cout << "Inhibitory columns must be less than SIZE, try again:" <<endl;
        cin >> inhibCols;
    }
    exCols = SIZE - inhibCols;
    inhibSparce = 0.4;
    exSparce = 0.1;
    gamma = 3.00;   // From Biology
    exConst = 1.054;  // What Guillaume sets it to in his paper
    inhibConst = -((SIZE-inhibCols)*gamma*exSparce*exConst)/(inhibSparce*inhibCols);
    //inhibConst = -1.9365;   //Valid for 200x200 cols, from Guillaume's paper
    
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
