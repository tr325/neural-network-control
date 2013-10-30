


/* Solves lyapunov eqn for 0 = L(P,A,U,s)   */ 
/* Returns pointer to P                     */
double* lyap_P(double *A[], int sizeA, double *U[], int sizeU);


/* Returns the smoothed spectral abcissa for W, when given a specific  */ 
/* epsilon to work towards                          */  
double simpleSSA(double *W[], double eps);
