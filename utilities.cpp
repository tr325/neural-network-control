/* File containing utility functions for Neural Network simulations   */

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include"utilities.h"

using namespace std; 


extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);
extern "C" double dgees_(char* jobvs, char* sort, bool* select, int* n, double* A, int* lda, int* sdim, double* wr, double* wi, double* vs, int* ldvs, double* work, int* lwork, bool* bwork, int* info);
extern "C" double dgemm_(char* trana, char* tranb, int* m, int* n, int* k, double* alpha, double* A, int* lda, double* B, int* ldb, double* beta, double* C, int* ldc);


/*  Takes inputs the required epsilon (larger->smoother SSA) and the input matrix A 
*   Calculates the value of the SSA using the lyapunov eqn solver  functions        */ 
double SimpleSSA(double *W[], double eps)
{
    double precision; 
    precision = 0.000001;
    
    // Generate A = (W-sI) 
    double *I[SIZE], *A[SIZE], *A_t[SIZE], *P[SIZE], *Q[SIZE];
    double *fP, *fQ, *fA, *fA_t;
    double spectralAbcissa, s, sOld;
    int l, loopcount;
    bool conv;
    
    conv = false;
    l = 0;
    loopcount = 0;
    s = 0.00;
    sOld = 0.01;
          
    for(int i=0; i<SIZE; i++)
    {
        I[i] = new double[SIZE];
        A[i] = new double[SIZE];
        A_t[i] = new double[SIZE]; 
        P[i] = new double[SIZE];
        Q[i] = new double[SIZE];
        I[i][i] = 1;                // The rest of the array default initialises to 0
    }
    
    
    // main loop for N-R root finding method
    while(!conv)
    {
        sOld = s; 
        
        for(int i=0; i<SIZE; i++)
        {
            for(int j=0; j<SIZE; j++)
            {
                A[i][j] = W[i][j] - s*I[i][j];
            }
        }
        // Form transpose(A)
        for(int i=0; i<SIZE; i++)
        {
            for(int j=0; j<SIZE; j++)
            {
                A_t[i][j] = A[j][i];
            }
        }     
                    
        // convert A and A_t to Schur form using dgees_ 
        fA = FArrayConvert(A, SIZE);
        Schur(&fA);     
        CArrayConvert(fA, A, SIZE); 
        
        fA_t = FArrayConvert(A_t, SIZE); 
        Schur(&fA_t); 
        CArrayConvert(fA_t, A_t, SIZE); 
        
        // Find SA value (max value of diagonals of Schur form) 
        // operates on first loop only
        if(l == 0)
        {
            spectralAbcissa = MaxDiag(A, SIZE);
            cout << setprecision(10) <<"Spectral abcissa value = " << spectralAbcissa <<endl; 
            l++;
        }            
        
        // Solve lyapunov equations for P and Q    
        fP = Lyap(fA, fA_t); 
        CArrayConvert(fP, P, SIZE); 
      
        fQ = Lyap(fA_t, fA); 
        CArrayConvert(fQ, Q, SIZE); 
       
        // use Newton-Raphson method to find new s value 
        conv = Newton(Q, P, eps, spectralAbcissa, precision, s);  
        
        //cout << "s = " << s << ", sOld = " << sOld <<endl;
        loopcount++;
       /* if(s == sOld)
        {
            cout << "desired SSA is below threshold limit" << endl;
            break;
        }*/
    }   
    
    cout << endl << "Smoothed spectral abcissa value = " << s << endl; 
    cout << "Spectral abcissa value = " << spectralAbcissa << endl; 
    cout << "Loop ran " << loopcount << " times" <<endl; 
    return s; 
}


/*  Finds the value of s for the next iteration of the loop     */
bool Newton(double *Q[], double *P[], double eps, double sA, double prec, double &s)
{
    double val, grad, e; 
    double *fP, *fQ, *fPQ;
    double *PQ[SIZE]; 
    
    for(int i=0; i<SIZE; i++)
    {
        PQ[i] = new double[SIZE];
    }    
    
    fP = FArrayConvert(P, SIZE); 
    fQ = FArrayConvert(Q, SIZE); 
    fPQ = FArrayConvert(PQ, SIZE);
    MatMult(fP, fQ, fPQ); 
    CArrayConvert(fPQ, PQ, SIZE); 
    
    val = Trace(Q, SIZE); 
    e = 1/eps;
  //  cout << "e = " << e << ", and val = " << val <<endl; 
    if(!(abs(val - e) < prec))
    {
        grad = -2*Trace(PQ, SIZE); 
        s = s + (1/grad)*(e - val);         
       // cout << (1/grad)*(e - val) << endl;
        if(s < sA)
        {
            s = sA + 0.1;
        }
       // cout << "gradient = " <<grad <<", s = " << s <<", val = " <<val << endl;  
        return false;
    }
    else 
    {
        return true; 
    } 
}    

/*  Performs matrix multiplication of matrices A and B,     
 *  storing the result in C                                 */
void MatMult(double A[], double B[], double C[])
{
    char TRANA, TRANB; 
    int M, N, K, LDA, LDB, LDC; 
    double ALPHA, BETA; 
    
    TRANA = 'N';
    TRANB = 'N';
    M = SIZE;
    N = SIZE;
    K = SIZE;
    LDA = SIZE;
    LDB = SIZE;
    LDC = SIZE;
    ALPHA = 1.00;
    BETA = 0.00;
    
    dgemm_(&TRANA, &TRANB, &M, &N, &K, &ALPHA, A, &LDA, B, &LDB, &BETA, C, &LDC);
    
    return;
}
    
    
/*  Calculates the trace of a matrix                    */
double Trace(double *A[], int dimA)
{
    double t = 0.00; 
    for(int i=0; i<dimA; i++)
    {
        t += A[i][i];
    }
    return t;
}

/*  Finds the maximum value of the diagonals of a matrix    */
double MaxDiag(double *A[], int dimA)
{
    double p = 0.00;
    for(int i=0; i<dimA; i++)
    {
        if(A[i][i] > p)
        {
            p = A[i][i];
        }
    }
    
    return p; 
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
            fArray[dimA*i + j] = A[j][i];
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

/* Populates a 2D array with the data held in an ascii file */ 
void loadMat(ifstream &file, double *A[], int matDim)
{
    if(file.is_open())
    {
        cout << "File is open" << endl;
        for(int i=0; i<matDim; i++)
        {
            for(int j=0; j<matDim; j++)
            {
                if(file.eof())
                {
                    cout << "ERROR: file ended before expected" <<endl; 
                    break;
                }
                file >> A[i][j];
            }
        } 
    }
    else
    {
        cout << "The file isn't open for reading" << endl;
    }
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
    
    //cout << "dtrsyl return status: " <<INFO << endl; 
    
    return fC;
}




