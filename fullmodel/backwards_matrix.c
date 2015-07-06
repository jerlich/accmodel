/*=================================================================
 * backwards_matrix.c
 *
 * Computes a one timestep backwards propagation matrix for Fokker-
 * Planck simulations. The equivalent Matlab code is
 *
 *      WB = W'.*(Pf(:,k-1)*ones(1,size(W,2)))./(ones(size(W,1),1)*Pf(:,k)');
 *
 * where W is the forward matrix and Pf(:,k) is the forward probability
 * distribution at timestep k. 
 *
 * Call as
 *
 *  >>  WB = backwards_matrix(W, Pfkm1, Pfk)
 *
 *
 *  To compile backwards_matrix.c, simply run, from within Matlab,
 *
 *     >> mex backwards_matrix.c

 *=================================================================*/

/* written by CDB March 2011 */
#include "mex.h"
#include "matrix.h"

void
        mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   int i, j;
   int M;
   int N;
   double *Bmatrix;
   double *Fmatrix;
   double *Fdistribm1;
   double *Fdistrib;
   
   /* Check for proper number of input and output arguments 
   if (nrhs != 3) {
      mexErrMsgTxt("Need three input args: W (forward matrix), "
              "Fkm1 (forward distrib at step k-1), Fk (forward distrib at step k).");
   }
   if(nlhs > 1){
      mexErrMsgTxt("Too many output arguments.");
   } */

   M = mxGetM(prhs[0]); N = mxGetN(prhs[0]);
   /* if ( M != N ) {
      mexErrMsgTxt("First arg must be a square matrix.");
   }
   if ( (mxGetM(prhs[1])!=M) || (mxGetN(prhs[1])!=1) ) {
      mexErrMsgTxt("Second arg must be a column vector, length same as columns of first arg");
   }
   if ( (mxGetM(prhs[2])!=M) || (mxGetN(prhs[2])!=1) ) {
      mexErrMsgTxt("Third arg must be a column vector, length same as columns of first arg");
   } */

   
   plhs[0]= mxCreateDoubleMatrix(M, N, mxREAL);
   
   Bmatrix = mxGetData(plhs[0]);
   
   Fmatrix    = mxGetData(prhs[0]);
   Fdistribm1 = mxGetData(prhs[1]);
   Fdistrib   = mxGetData(prhs[2]);

   for ( i=0; i<M; i++ ) {
      for ( j=0; j<N; j++ ) {
         Bmatrix[j*M+i] = Fmatrix[i*M+j]*Fdistribm1[i]/Fdistrib[j];
      }
   }
   
      
}

