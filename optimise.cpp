//  File to test optimization of W matrix, using SSA solving methods in utilities.cpp


#include<iostream>
#include<cmath>
#include<cstdlib>
#include"optimise.h"
#include"utilities.h"
#include<vector>
#include<ctime>

using namespace std;


/*  Takes the input weighting matrix, and optimises the second half
 *  to minimise the ssa value.  Optimisation is by gradient descent */
void OptimiseWMat(double *W[], int *B[], double eps, int inhibNum, bool wPlast, int SIZE)
{
    double      *A[SIZE];
    double      *I[SIZE];
    double      *P[SIZE];
    double      *Q[SIZE];
    double      *V[SIZE];
    double      *S_Vect[SIZE];
    double      *S_Vect_T[SIZE];
    double      *gradMat[SIZE];
    double      *QP[SIZE];
    double      *origW[SIZE];
    bool        conv;
    double      *sa;
    double      ssa;
    double      ssaOld;
    double      precision;
    double      descentRate;
    double      avgStrengh;
    int         loopcount;
    int         numReformed;
    int         itLimit;
    vector<double>  ssaWindow;
    clock_t     startTime;

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
        origW[i] = new double[SIZE];
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
            origW[i][j] = W[i][j];
        }
    }

    conv = false;
    ssa = 0.0;
    ssaOld = 100.0;  //initialised to greater than ssa for descentRate manipulation
    precision = 0.001;
    loopcount = 0;
    descentRate = 10;   // determines how far down slope each iteration updates
                        // is varied when unstable behaviour occurs
    itLimit = 400;
    Reparam(W, B, V, SIZE);
    //~ OutputMat("reparamW.ascii", W, SIZE);
    
    ofstream resFile;       //stores ssa values over loop
    ofstream recFile;         // stores ssa start and end values
    ofstream reformRecord;  // stores record of reform events
    resFile.open("SsaOptimise.ascii");
    reformRecord.open("refSynRecord.ascii");
    recFile.open("SSArecord.ascii");

    
    // loops until the gradient descent algorithm has converged
    while(!conv)
    {
        loopcount++;
        
        //cout << "Optimisation loopcount = " <<loopcount << endl; 
        // solve SSA for current W matrix
        ssaOld = ssa;
        ssa = SimpleSSA(W, P, Q, S_Vect, eps, SIZE);

        // window for testing convergence
        ssaWindow.push_back(ssa);
        if(ssaWindow.size() > 100)
        {
            ssaWindow.erase(ssaWindow.begin());
            if(abs(ssaWindow.front() - ssaWindow.back()) < precision)
            {
                cout <<"Window difference measure = " <<abs(ssaWindow.front() - ssaWindow.back()) <<endl;
                cout << "Converged on loop " <<loopcount <<endl;
                conv = true;
            }
        }
        
        //~ if(loopcount == itLimit)
        //~ {
            //~ cout <<"Iteration limit" <<endl;
            //~ conv = true;
        //~ }
        
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
        if(wPlast)
        {
            numReformed = ReformSyn(V, B, inhibNum, SIZE);
        }
        RecalcW(W, B, V, SIZE);

        
        if(numReformed)
        {
            reformRecord << loopcount <<"   " <<numReformed <<endl;
        }            
        
        if(loopcount == 1)
        {
            recFile << ssa << "   ";
        }
        
        int d = loopcount%200;
        if(d == 1)
        {
            cout <<"SSA on loop " << loopcount << ": " <<ssa <<endl;
        } 
        
        resFile << ssa <<endl;

        if(ssa > 100)
        {
            cout <<"Unstable behaviour" <<endl;
            for(int i=0; i<SIZE; i++)
            {
                for(int j=0; j<SIZE; j++)
                {
                    W[i][j] = origW[i][j];
                }
            }
            loopcount = 0;
            descentRate = descentRate/1.5;
            itLimit = int(ceil(itLimit*1.5));
        }

        //cout <<endl <<endl;
        //~ cout <<"SSA on loop " << loopcount << ": " <<ssa <<endl; 
        
    }
    // records final SSA for octave plotting
    cout <<"Final SSA = " <<ssa <<endl;
    recFile <<ssa <<endl;
    recFile.close();

    resFile.close();
    reformRecord.close();

    ofstream loopCountRecord;
    loopCountRecord.open("OptimisationCovergenceTime.ascii");
    loopCountRecord << "SIZE = " <<SIZE << "no. of inhibCols = " <<inhibNum <<", no. of loops to stabilise = " <<loopcount << endl;
    loopCountRecord.close();
    //~ OutputMat("reparStabW.ascii", W, SIZE);
    //~ OutputMat("reparStabV.ascii", V, SIZE);
    //~ OutputMat("reparStabB.ascii", B, SIZE);
    
    for(int i=0; i<SIZE; i++)
    {
        delete[]  A[i];
        delete[]  I[i];
        delete[]  P[i];
        delete[]  Q[i];
        delete[]  V[i];
        delete[]  S_Vect[i];
        delete[]  S_Vect_T[i];
        delete[]  gradMat[i];
        delete[]  QP[i];
        delete[]  origW[i];
    }

    return;
}


/*  Forms new synapses if some decay to zero        */
int ReformSyn(double *V[], int *B[], int inhibNum, int SIZE)
{
    double initVal;
    double decayVal;
    int exNum;
    int newCol;
    int fInd;
    int refNum;
    bool test;
    vector<int> freeCols;
    
    exNum = SIZE - inhibNum;
    initVal = -3;
    decayVal = -3;
    refNum = 0;
    
    //~ OutputMat("refStabTestV1.ascii", V, SIZE);
    //~ OutputMat("refStabTestB1.ascii", B, SIZE);
    
    srand(time(NULL));  //seeds the rand() function with the UNIX time

    
    //cout << "reforming..." <<endl;
    for(int i=0; i<SIZE; i++)
    {
        for(int j=exNum; j<SIZE; j++)
        {
            if(V[i][j] < decayVal)
            {
                //~ cout << "Decayed synapse: " <<i+1 <<", " <<j+1 <<"; strength = "<<V[i][j] <<endl;
                for(int k=exNum; k<SIZE; k++)
                {
                    if((B[i][k] == 0) && (i != k))
                    {
                        freeCols.push_back(k);
                    }
                }

                fInd = rand()%freeCols.size();
                newCol = freeCols[fInd];
                //~ cout << "Formed synapse: " <<i+1 <<", " <<newCol+1 <<"; value updated from " <<V[i][newCol] <<" to ";
                B[i][newCol] = -1;
                V[i][newCol] = initVal;
                //~ cout <<V[i][newCol] <<endl;
                V[i][j] = 0;
                B[i][j] = 0;
                refNum ++;  
                test = true;           



                //~ newCol = (rand() % inhibNum) + exNum;
                //~ while((B[i][newCol] != 0) || (newCol == i))
                //~ {
                    //~ newCol = (rand() % inhibNum) + exNum;
                //~ }
                
            }
        }
    } 
    
    
    //~ if(test)
    //~ {
        //~ OutputMat("refStabTestV2.ascii", V, SIZE);
        //~ OutputMat("refStabTestB2.ascii", B, SIZE);
    //~ }
    
    return refNum;
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

