/*  Calculates the gradient matrix by numerical methods */

#include<iostream>
#include"utilities.h"


using namespace std;

int main()
{

    int** m = new int*[2];
    m[0] = new int[2];
    m[1] = new int[2];
    
    m[0][0] = 1; 
    m[0][1] = 2; 
    m[1][0] = 3;
    m[1][1] = 12;
    
    double t; 
    
    t = Trace(m, 2);
    
    cout <<t <<endl;
    
    return 0;
}
