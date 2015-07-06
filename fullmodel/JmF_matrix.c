/*=================================================================
 * JmF_matrix.c
 *
 * Computes the J_ij/F_ij matrix
 *
 *      JmF = (ones(n,1)*Pf(:,k-1)')./(Pf(:,k)*ones(1, n)).*(Pb(:,k)*ones(1, n));
 *
 * where Pf(:,k) is the forward probability distribution at timestep k and 
 * Pb(:,k) is the backward probability distribution at timestep k
 * 
 *
 * Call as
 *
 *  >>  JmF = JmF_matrix(Pfkm1, Pfk, Pbk)
 *
 *
 *  To compile backwards_matrix.c, simply run, from within Matlab,
 *
 *     >> mex JmF_matrix.c

 *=================================================================*/

/* written by BWB August 2011 */
#include "mex.h"
#include "matrix.h"

void
        mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   int i, j;
   int M;
   double *JmFmatrix;
   double *Pfkm1;
   double *Pfk;
   double *Pbk;
   
   M = mxGetM(prhs[0]);
   
   plhs[0]= mxCreateDoubleMatrix(M, M, mxREAL);
   
   JmFmatrix = mxGetData(plhs[0]);
   
   Pfkm1 = mxGetData(prhs[0]);
   Pfk   = mxGetData(prhs[1]);
   Pbk   = mxGetData(prhs[2]);

   for ( i=0; i<M; i++ ) {
      for ( j=0; j<M; j++ ) {
         JmFmatrix[j*M+i] = Pfkm1[j]/Pfk[i]*Pbk[i];
      }
   }
   
      
}

