/*  tests the MatMult function */
/*  tests the MatMult function */

#include<fstream>
#include<iostream>
#include"utilities.h"
#include"optimise.h"


using namespace std;

int main()
{
    ifstream pfile;
    ifstream qfile;
    ifstream qpfile;
    ofstream opfile;
    int SIZE;

    pfile.open("TESTP.ascii");
    qfile.open("TESTQ.ascii");
    qpfile.open("TESTQP.ascii");
    SIZE = MatSize(pfile);
    
    double *mat1[SIZE];
    double *mat2[SIZE];
    double *ref[SIZE];
    double *ans[SIZE];
    double *f1;
    double *f2;
    double *fR;
    double *fA;
    
    for(int i=0; i<SIZE; i++)
    {
        mat1[i] = new double[SIZE];
        mat2[i] = new double[SIZE];
        ans[i] = new double[SIZE];
        ref[i] = new double[SIZE];
    }
    
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            qfile >> mat1[i][j];
            pfile >> mat2[i][j];
            qpfile >> ans[i][j];
            ref[i][j] = mat1[i][j];
        }
    }
        
    fR = FArrayConvert(ref, SIZE);
    CArrayConvert(fR, ref, SIZE);
    
    double diff;
    //cout << "converted there and back:" <<endl;
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            diff = ref[i][j] - mat1[i][j];
            if(diff != 0)
            {
                cout << "out by " << diff <<endl;
            }
        }
    }
        
        
        
        
        
        
    f1 = FArrayConvert(mat1, SIZE);
    f2 = FArrayConvert(mat2, SIZE);
    fA = FArrayConvert(ans, SIZE);
    MatMult(f1, f2, fA, SIZE);
    
    CArrayConvert(fA, ans, SIZE);

    cout << "fA:" <<endl;
    
    for(int i=0; i<SIZE*SIZE; i++)
    {
        cout << fA[i] <<"  ";
        
    }
    
    cout <<endl <<endl;

    cout << "Result" <<endl;
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            cout << ans[i][j] <<"  ";
        }
        cout <<endl;
    }    
    
    
    opfile.open("TESTQP.ascii");
    for(int i=0; i<SIZE; i++)
    {
        for(int j=0; j<SIZE; j++)
        {
            opfile << ans[i][j] <<"  ";
        }
        opfile <<endl;
    }    

    return 0;
    
}    
