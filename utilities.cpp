/* File containing utility functions for Neural Network simulations   */

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<cctype>
#include"utilities.h"

using namespace std; 


extern "C" double dtrsyl_(char* trana, char* tranb, int* isgn, int* orda, int* ordb, double* A, int* lda, double* B, int* ldb, double* C, int* ldc, double* scale, int* info);
extern "C" double dgees_(char* jobvs, char* sort, bool* select, int* n, double* A, int* lda, int* sdim, double* wr, double* wi, double* vs, int* ldvs, double* work, int* lwork, bool* bwork, int* info);
extern "C" double dgemm_(char* trana, char* tranb, int* m, int* n, int* k, double* alpha, double* A, int* lda, double* B, int* ldb, double* beta, double* C, int* ldc);


/*  Takes inputs the required epsilon (larger->smoother SSA) and the input matrix A 
*   Calculates the value of the SSA using the lyapunov eqn solver functions         */ 
double SimpleSSA(double *W[], double *P[], double *Q[], double eps, int SIZE)
{
    double precision; 
    precision = 0.000001;
    
    double *I[SIZE], *A[SIZE];
    double *fP, *fQ, *fA;
    double spectralAbcissa, s;
    int loopcount;
    bool conv;
    
    conv = false;
    loopcount = 0;
    s = 0;
          
    for(int i=0; i<SIZE; i++)
    {
        I[i] = new double[SIZE];
        A[i] = new double[SIZE];
        I[i][i] = 1;     // The rest of the array default initialises to 0
    }
    
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            A[i][j] = W[i][j];
        }
    }
    
        
    // main loop for N-R root finding method
    while(!conv)
    {    
        // Generate A = (W-sI)     
        for(int i=0; i<SIZE; i++)
        {
            for(int j=0; j<SIZE; j++)
            {
                A[i][j] = W[i][j] - s*I[i][j];
            }
        }
        
        // convert A and A_t to Schur form using dgees_
        fA = FArrayConvert(A, SIZE);
        Schur(fA, SIZE);
        CArrayConvert(fA, A, SIZE);

        // Find SA value (max value of diagonals of Schur form) 
        // operates on first loop only
        if(loopcount == 0)
        {
            spectralAbcissa = MaxDiag(A, SIZE);
            //cout << setprecision(10) <<"Spectral abcissa value = " << spectralAbcissa <<endl; 
            s = spectralAbcissa + 0.1;
            //cout << "Initial s =" << s << endl;
            /*
            ofstream opfile;
            opfile.open("TESTPrintA");
            for(int i=0; i<SIZE; i++)
            {
                for(int j=0; j<SIZE; j++)
                {
                    opfile << A[i][j] << " ";
                }
                opfile <<endl;
            }
            * */
        }           
        else
        {
            // Solve lyapunov equations for P and Q            
            fQ = Lyap(fA, true, SIZE);
            CArrayConvert(fQ, Q, SIZE);
            delete[] fQ;
            fP = Lyap(fA, false, SIZE);
            CArrayConvert(fP, P, SIZE);
            delete[] fP;
            
            // use Newton-Raphson method to find new s value            
            conv = Newton(Q, P, eps, spectralAbcissa, precision, s, SIZE);
        }
        
        delete[] fA;
        loopcount++;
    }
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] I[i];
        delete[] A[i];
    }
    
    cout << "Smoothed spectral abcissa value = " << s << endl;
    cout << "Spectral abcissa value = " << spectralAbcissa << endl;
    cout << "Loop ran " << loopcount << " times" <<endl <<endl;
    return s;
}

/*  Finds the value of s for the next iteration of the loop     */
bool Newton(double *Q[], double *P[], double eps, double sA, double prec, double &s, int SIZE)
{
    double val, grad, e;
    double *fP, *fQ, *fPQ;
    double *PQ[SIZE];
    bool conv; 
    
    for(int i=0; i<SIZE; i++)
    {
        PQ[i] = new double[SIZE];
    }
    
    val = Trace(Q, SIZE);
    
    fP = FArrayConvert(P, SIZE);
    fQ = FArrayConvert(Q, SIZE);
    fPQ = FArrayConvert(PQ, SIZE);
    MatMult(fP, fQ, fPQ, SIZE);
    CArrayConvert(fPQ, PQ, SIZE);
    delete[] fP;
    delete[] fPQ;
    delete[] fQ;

    e = 1/eps;
    //cout << "e = " << e << ", and val = " << val <<endl; 
    if(!(abs(val - e) < prec))
    {
        grad = -2*Trace(PQ, SIZE);
        s = s + (1/grad)*(e - val);
        
        if(s < sA + 0.00001)
        {
            s = sA + 0.00001;
        }
        //cout << "gradient = " <<grad <<", s = " << s <<", val = " <<val << endl;
        //cout << "s jump " << (1/grad)*(e - val) << endl;
        conv = false;
    }
    else
    {
        conv = true;
    }
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] PQ[i];
    }
    
    return conv; 
}    

/*  Performs matrix multiplication of matrices A and B,     
 *  storing the result in C                                 */
void MatMult(double A[], double B[], double C[], int SIZE)
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

/*  Converts input matrix A to Schur form (req for dtrsyl)  
 *  input matrix should be coverted for passing to Fortran by 
 *  FArrayConvert() before being passed to this funciton        */
void Schur(double A[], int SIZE)
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
    
    dgees_(&JOBVS, &SORT, 0, &N, A, &LDA, &SDIM, *EIG_Re, *EIG_Im, *VS, &LDVS, *WORK, &LWORK, *BWORK, &INFO);
    
    return;    
}

/* Populates a 2D array with the data held in an ascii file */ 
void loadMat(ifstream &file, double *A[], int SIZE)
{
    int j=0; 
    int i=0;

    if(file.is_open())
    {
        
        //cout << "File is open" << endl;

        for(i=0; i<SIZE; i++)
        {
            for(j=0; j<SIZE; j++)
            {
                if(file.eof())
                {
                    cout << "ERROR: file ended before expected" <<endl; 
                    break;
                }
                //cout << "pleep" <<endl; 
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

/* Returns dimension of input matrix from the file  */
int MatSize(ifstream &file)
{
    int numlines = 0;
    char c;
    string s;
    while(!file.eof())
    {
        getline(file, s);
        numlines++;
    }
    
    numlines--;  // the eof char counts as an extra line
    
    file.clear();
    file.seekg(0);

    return numlines;
}



/* Solves lyapunov eqn for 0 = L(P,A,U,s)       */ 
/* Returns pointer to P                         */
/* Solves: (W-sI)P + P(W-sI)^T = -I             */
/* Form:  A*P + P*B = C                         */
double* Lyap(double A[], bool TRAN, int SIZE)
{
    char TRANA, TRANB;
    int ISGN, ORDA, ORDB, LDA, LDB, LDC, INFO;
    
    if(TRAN)
    {
        TRANA = 'T';
        TRANB = 'N';
    }
    else
    {
        TRANA = 'N';
        TRANB = 'T';
    }
    ISGN = 1;

    ORDA = SIZE;
    //ORDB = SIZE;
    LDA = SIZE;
    //LDB = SIZE;
    LDC = SIZE;
    
    double *C[SIZE];
    double *fC;
    
    for(int i=0; i<SIZE; i++)
    {
        C[i] = new double[SIZE];
        C[i][i] = -1.00;
    }

    fC = FArrayConvert(C, SIZE);

    double SCALE = 1.00;
    
    dtrsyl_(&TRANA, &TRANB, &ISGN, &ORDA, &ORDA, A, &LDA, A, &LDA, fC, &LDC, &SCALE, &INFO);  
    
    for(int i=0; i<SIZE*SIZE; i++)
    {
        fC[i] = fC[i]/SCALE;
    }
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] C[i];
    }
    
    //cout << "dtrsyl return status: " <<INFO << endl; 
    
    return fC;
}
