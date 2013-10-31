/* File containing utility functions for Neural Network simulations   */

#include<iostream>
#include"utilities.h"

using namespace std; 



extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);

extern "C" double dgees_(char* jobvs, char* sort, bool* select, int* n, double* A, int* lda, int* sdim, double* wr, double* wi, double* vs, int* ldvs, double* work, int* lwork, bool* bwork, int* info);




/*  Takes inputs the required epsilon (larger->smoother SSA) and the input matrix A  */ 
/*  Calculates the value of the SSA using the lyapunov eqn solver  functions  */ 
double SimpleSSA(double *W[], double eps)
{
    
    // Generate A = (W-sI) 
    double *I[SIZE];
    double *A[SIZE];
    double *A_t[SIZE]; 
      
    for(int i=0; i<SIZE; i++)
    {
        I[i] = new double[SIZE];
        A[i] = new double[SIZE];
        A_t[i]= new double[SIZE]; 
        I[i][i] = 1;                // The rest of the array default initialises to 0
        for(int j=0; j<SIZE; j++)
        {
            A[i][j] = W[i][j] - eps*I[i][j];
        }
    }
    
    cout << "Matrix A:" <<endl; 
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            cout << A[i][j] << "    ";
        }
        cout << endl; 
    }
    cout << endl;
    
    // Form transpose(A)
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            A_t[i][j] = A[j][i];
        }
    }     
                
    // convert A and A_t to Schur form using dgees_ 
    double *fA;
    double *fA_t; 
    
    fA = FArrayConvert(A, SIZE);
    Schur(&fA);     
    CArrayConvert(fA, A, SIZE); 
    
    fA_t = FArrayConvert(A_t, SIZE); 
    Schur(&fA_t); 
    CArrayConvert(fA_t, A_t, SIZE); 
    
    cout << endl << "Schur(A): " << endl;
     for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            cout << A[i][j] << "    ";
        }
        cout << endl; 
    }
    cout << endl;
    
    // Solve lyapunov equations for P and Q 
    double *fP, *fQ;
    double *P[SIZE], *Q[SIZE];
    for(int i=0; i<SIZE; i++)
    {
        P[i] = new double[SIZE];
        Q[i] = new double[SIZE];
    }
     
    fP = Lyap(fA, fA_t);
 //   for(int i=0; i<(SIZE*SIZE); i++)
 //       cout << fP[i] << endl; 
    CArrayConvert(fP, P, SIZE); 
    
    cout << "Matrix P:" <<endl; 
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            cout << P[i][j] << "    ";
        }
        cout << endl; 
    }
    cout << endl;
    
    fQ = Lyap(fA_t, fA); 
    CArrayConvert(fQ, Q, SIZE); 
    
    cout << "Matrix Q:" <<endl; 
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            cout << Q[i][j] << "    ";
        }
        cout << endl; 
    }
    cout << endl;
    
    
    // use Newton-Raphson method to locate SSA   
    
}

/*  Converts input array to Column-major 1-D array for      */
/*  passing to Fortran subroutines                          */
double* FArrayConvert(double *A[], int dimA)
{
    double *fArray;
    fArray = new double[dimA*dimA]; 
    for(int i=0; i<dimA; i++)
    {
        for(int j=0; j<dimA; j++)
        {
            fArray[3*i + j] = A[j][i];
        }
    }
    
    return fArray;    
}

/*  Repopulates array cA with the results of a Fortran  */
/*  subroutine, held in array fA                        */
void CArrayConvert(double fA[], double *cA[], int dimCA)
{
    for(int i=0; i<dimCA; i++)
    {
        for(int j=0; j<dimCA; j++)
        {
            cA[j][i] = fA[dimCA*i + j];
        }
    }
}

/*  Converts input matrix A to Schur form (req for dtrsyl)  */
void Schur(double *A[])
{
    char JOBVS, SORT;
    int N, LDA, LDVS, LWORK, SDIM, INFO;  
    double *EIG_Re[SIZE], *EIG_Im[SIZE], *VS[SIZE], *WORK[SIZE]; 
    bool *BWORK[SIZE]; 
    
    for(int i=0; i<SIZE; i++)
    {
        EIG_Re[i] = new double[SIZE];
        EIG_Im[i] = new double[SIZE];
        VS[i] = new double[SIZE];
        WORK[i] = new double[SIZE];
        BWORK[i] = new bool[SIZE];
    }
    JOBVS = 'N';
    SORT = 'N';
    N = SIZE;
    LDA = SIZE; 
    LDVS = SIZE;
    LWORK = 3*SIZE;
    
    dgees_(&JOBVS, &SORT, 0, &N, *A, &LDA, &SDIM, *EIG_Re, *EIG_Im, *VS, &LDVS, *WORK, &LWORK, *BWORK, &INFO);
    
    //cout << "dgees_ return status: " << INFO <<endl; 
    
    return;    
}


/* Solves lyapunov eqn for 0 = L(P,A,U,s)   */ 
/* Returns pointer to P                     */
/* Solves: (W-sI)P + P(W-sI)^T = -2I        */
/* Form:  A*P + P*B = C                     */
double* Lyap(double A[], double B[])
{
    char TRANA, TRANB;
    int ISGN, ORDA, ORDB, LDA, LDB, LDC, INFO;
    
    TRANA = 'N';
    TRANB = 'N';
    ISGN = 1;

    ORDA = SIZE;
    ORDB = SIZE;
    LDA = SIZE;
    LDB = SIZE;
    LDC = SIZE;
    
    double *C[SIZE];
    double *fC;
    for(int i=0; i<SIZE; i++)
    {
        C[i] = new double[SIZE];
        C[i][i] = -2.00;
    }
            
    fC = FArrayConvert(C, SIZE);
    
    double SCALE = 1.00;
    
    dtrsyl_(&TRANA, &TRANB, &ISGN, &ORDA, &ORDB, A, &LDA, B, &LDB, fC, &LDC, &SCALE, &INFO);  
    
    cout << "dtrsyl return status: " <<INFO << endl; 
    
    return fC;
}




