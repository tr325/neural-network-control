//  File to test optimization of W matrix, using SSA solving methods in utilities.cpp


#include<iostream>
#include<cmath>
#include<cstdlib>
#include"optimise.h"
#include"utilities.h"

using namespace std;


/*  Takes the input weighting matrix, and optimises the second half
 *  to minimise the ssa value.  Optimisation is by gradient descent */
void OptimiseWMat(double *W[], int *B[], double eps, int inhibNum, int SIZE)
{
    double  *A[SIZE];
    double  *I[SIZE];
    double  *P[SIZE];
    double  *Q[SIZE];
    double  *V[SIZE];
    double  *S_Vect[SIZE];
    double  *S_Vect_T[SIZE];
    double  *gradMat[SIZE];
    double  *QP[SIZE];
    bool    conv;
    double  *sa;
    double  ssa;
    double  ssaOld;
    double  precision;
    double  descentRate;
    double  avgStrengh;
    int     loopcount;

    for(int i=0; i<SIZE; i++)
    {
        A[i] = new double[SIZE];
        I[i] = new double[SIZE];
        P[i] = new double[SIZE];
        Q[i] = new double[SIZE];
        V[i] = new double[SIZE];
        S_Vect[i] = new double[SIZE];
        S_Vect_T[i] = new double[SIZE];
        gradMat[i] = new double[SIZE];
        QP[i] = new double[SIZE];
        for(int j=0; j<SIZE; j++)
        {
            if(i==j)
            {
                I[i][j] = 1;
            }
            else
            {
                I[i][j] = 0;
            }
        }
    }

    conv = false;
    ssa = 0.0;
    ssaOld = 10.0;  //initialised to greater than ssa for descentRate manipulation
    precision = 0.0000001;
    loopcount = 0;
    descentRate = 5;   // determines how far down slope each iteration updates
    Reparam(W, B, V, SIZE);
    OutputMat("reparamW.ascii", W, SIZE);
    ofstream resFile;
    resFile.open("SsaOptimise.ascii");
    
    
    // loops until the gradient descent algorithm has converged
    while(!conv)
    {
        //cout << "Optimisation loopcount = " <<loopcount << endl; 
        // solve SSA for current W matrix
        ssaOld = ssa;
        ssa = SimpleSSA(W, P, Q, S_Vect, eps, SIZE);
        
        for(int i=0; i<SIZE; i++)
        {
            for(int j=0; j<SIZE; j++)
            {
                A[i][j] = W[i][j] - ssa*I[i][j];
                S_Vect_T[i][j] = S_Vect[j][i];
            }
        }
        
        // find the gradient matrix, d(SSA)/dW;
        MatMult(Q, P, QP, SIZE);
        FormGradMat(gradMat, W, QP, S_Vect, S_Vect_T, ssa, SIZE);
        
        // perform gradient descent to optimise V, then recalculate W
        GradDescent(V, gradMat, descentRate, inhibNum, SIZE);
        ReformSyn(V, B, inhibNum, SIZE);
        RecalcW(W, B, V, SIZE);
        
        
        if(abs(ssa - ssaOld) < precision)
        {
            cout << "Converged" <<endl;
            conv = true;
        }
        /*
        if(loopcount == 1000)
        {
            cout <<"That's enough i think..." <<endl;
            conv = true;
        }
        */
        if(ssa < 1.00)
        {
            cout << "Stabilised" << endl;
            conv = true;
        }
        if(ssa > 100)
        {
            cout <<"Unstable behaviour" <<endl;
            OutputMat("divergentW.ascii", W, SIZE);
            OutputMat("divergentB.ascii", B, SIZE);
            OutputMat("divergentV.ascii", V, SIZE);
            exit(0);
        }
        
        resFile << ssa <<endl;
                
        //cout <<endl <<endl;
        cout <<"SSA on loop " << loopcount << ": " <<ssa <<endl; 
        loopcount++;
    }
    resFile.close();

    OutputMat("reparStabW.ascii", W, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        delete[] P[i];
        delete[] Q[i];
        delete[] gradMat[i];
        delete[] QP[i];
        delete[] A[i];
        delete[] I[i];
        delete[] V[i];
        delete[] S_Vect_T[i];
    }

    return;
}


/*  Forms new synapses if some decay to zero        */
void ReformSyn(double *V[], int *B[], int inhibNum, int SIZE)
{
    double initVal;
    double decayVal;
    int exNum;
    int newCol;
    bool test; 
    
    exNum = SIZE - inhibNum;
    initVal = -5;
    decayVal = -5;
    
    OutputMat("refStabTestV1.ascii", V, SIZE);
    OutputMat("refStabTestB1.ascii", B, SIZE);
    
    //cout << "reforming..." <<endl;
    for(int i=0; i<SIZE; i++)
    {
        for(int j=exNum; j<SIZE; j++)
        {
            if(V[i][j] < decayVal)
            {
                cout << "Decayed synapse: " <<i+1 <<", " <<j+1 <<"; strength = "<<V[i][j] <<endl;
                newCol = (rand() % inhibNum) + exNum;
                while((B[i][newCol] != 0) || (newCol == i))
                {
                    newCol = (rand() % inhibNum) + exNum;
                }
                cout << "Formed synapse: " <<i+1 <<", " <<newCol+1 <<"; value updated from " <<V[i][newCol] <<" to ";
                B[i][newCol] = -1;
                V[i][newCol] = initVal;
                cout <<V[i][newCol] <<endl;
                V[i][j] = 0;
                B[i][j] = 0;
                
                test = true;
            }
        }
    } 
    
    
    if(test)
    {
        OutputMat("refStabTestV2.ascii", V, SIZE);
        OutputMat("refStabTestB2.ascii", B, SIZE);
    }
    
    return;
}


/*  Forms the gradient matrix. 
 *  CURRENTLY optimising over V - so gradMat = d(ssa)/dV        */
void FormGradMat(double *gradMat[], double *W[], double *QP[], double *S_Vect[], double *S_Vect_T[], double ssa, int SIZE)
{
    double *G[SIZE];
    double *foo[SIZE];
    double tr;
    
    for(int i=0; i<SIZE; i++)
    {
        G[i] = new double[SIZE];
        foo[i] = new double[SIZE];
    }
    
    //cout << "Forming gradmat..." << endl;
    
    //  Calculate the matrix of schur vectors, V
    //~ Schur(A, V, SIZE);
    //~ for(int i=0; i<SIZE; i++)
    //~ {
        //~ for(int j=0; j<SIZE; j++)
        //~ {
            //~ Vt[i][j] = V[j][i];
        //~ }
    //~ }
    tr = Trace(QP, SIZE);    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            G[i][j] = QP[i][j]/tr;
        }
    }
    
    //  Normalise by multiplying by Vector matrix and its transpose(pre- and post-)
    MatMult(S_Vect, G, foo, SIZE);
    MatMult(foo, S_Vect_T, gradMat, SIZE);
    
    //  Adjust by chain rule, such that GradMat contains gradients wrt V, not W
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            // Forms gradmat for steepest gradient descent
            gradMat[i][j] *= W[i][j];
        }
    }

    for(int i=0; i<SIZE; i++)
    {
        delete[] G[i];
        delete[] foo[i];
    }
    
    //cout <<"exit of FormGradMat()" <<endl; 
    
    return;
}


/*  Performs gradient descent on V, based on gradients in gradMat  */
void GradDescent(double *V[], double *gradMat[], double descentRate, int inhibNum, int SIZE)
{
    int exNum; 
    exNum = SIZE - inhibNum; 
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=(exNum); j<SIZE; j++)
        {
            if(i != j)
            {
                V[i][j] -= descentRate*gradMat[i][j];
            }
        }
    }
    return;
}


/*  Recalculates the W matrix after gradient descent manipulation of V  */
void RecalcW(double *W[], int *B[], double *V[], int SIZE)
{
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            W[i][j] = B[i][j]*exp(V[i][j]);
        }
    }
    
    return;    
}


/*  From an input W matrix, finds the corresponding V matrix. 
 *  w_ij = b_ij*exp(v_ij)   ->    v_ij = ln(w_ij/b_ij)  if bij != 0 */ 
void Reparam(double *W[], int *B[], double *V[], int SIZE)
{
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            if(B[i][j] != 0)
            {
                //cout<<"nonzero Bij at "<<i<<"  " <<j<<endl;
                V[i][j] = log(W[i][j]/B[i][j]);
            }
            else
            {
                //cout << "i = " <<i << ", j = " <<j <<endl;
                V[i][j] = 0;
            }
        }
    }
    return;
}

