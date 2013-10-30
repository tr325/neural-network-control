/* File containing utility functions for Neural Network simulations   */


#include<iostream>
#include"utilities.h"

using namespace std; 

extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);


/*  Takes inputs the required epsilon (larger->smoother SSA) and the input matrix A  */ 
/*  Calculates the value of the SSA using the lyapunov eqn solver  functions  */ 
double simpleSSA(double* W[], double eps, int dimW)
{
    // Generate A = (W-sI) 
    double *I[dimW];
    for(int i=0; i<dimW; i++)
    {
        I[i] = new double[eps]();
        I[i][i] = 1;
    }  
   
    // convert to Schur form using dgees_ 
   
    // Solve lyapunov equations for P and Q 
   
    // use Newton-Raphson method to locate SSA   
}




/* Solves lyapunov eqn for 0 = L(P,A,U,s)   */ 
/* Returns pointer to P                     */
/* Solves: (W-sI)P + P(W-sI)^T = -2I        */
double* lyap_P(double *A[], int sizeA, double *B[], int sizeB)
{
    char TRANA, TRANB;
    int ISGN, ORDA, ORDB, LDA, LDB, LDC, INFO;
    
    TRANA = 'N';
    TRANB = 'N';
    ISGN = 1;

    ORDA = sizeA;
    ORDB = sizeB;
    LDA = sizeA;
    LDB = sizeB;
    LDC = sizeA;
    
    double *C[3];
    for(int i=0; i<3; i++)
    {
        C[i] = new double[3];
    }
    double SCALE = 2.00;
    
    dtrsyl_(&TRANA, &TRANB, &ISGN, &ORDA, &ORDB, *A, &LDA, *B, &LDB, *C, &LDC, &SCALE, &INFO);  
    
    /*
    \\ Look, I can pass and access 2d matrices! 
    double *mat[3];
    
    for(int i=0; i<3; i++)
    {
        mat[i] = new double[3];
    }
    for(int j=0; j<sizeA; j++)
    {	           
        for(int k=0; k<sizeA; k++)
        {    
            mat[j][k] = A[j][k];
        }
    }
    
    cout << mat[1][1] <<endl; 
*/    //return mat;
}




