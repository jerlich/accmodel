/*=================================================================
 * make_F_struct.c
 *
 *
 * Call as
 *
 *  >>  WB = make_F_struct(??)
 *
 *
 *  To compile backwards_matrix.c, simply run, from within Matlab,
 *
 *     >> mex make_F_struct.c

 *=================================================================*/

/* written by BWB Aug 2011 */
#include "mex.h"
#include "matrix.h"
//#define _USE_MATH_DEFINES // for C
#include <math.h>
#include <String.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   int i, j, k;
   int N;
   double *x, *dB, dx;
   double sigma2;
   double lambda;
   double h;
   double dt;
   
   mxArray *W;
   
   int ndeltas;
   char **fieldnames;

   double *F;
   double *dFds;
   double *dFdsig2;
   double *dFdB;
   double *dFdl;
   double *dFdh;
   mxArray *mxF;
   mxArray *mxdFds;
   mxArray *mxdFdsig2;
   mxArray *mxdFdB;
   mxArray *mxdFdl;
   mxArray *mxdFdh;
   
  
   /* inputs */
   N  = mxGetNumberOfElements(prhs[0]);
   x  = mxGetData(prhs[0]);
   dB = mxGetData(prhs[1]);
   sigma2 = mxGetScalar(prhs[2]);
   lambda = mxGetScalar(prhs[3]);
   h      = mxGetScalar(prhs[4]);
   dt     = mxGetScalar(prhs[5]);
   dx     = mxGetScalar(prhs[6]);   
   ///* ouputs */
   fieldnames[0]="F";
   fieldnames[1]="dFdsig2";
   fieldnames[2]="dFdB";
   fieldnames[3]="dFdl";
   fieldnames[4]="dFdh";
    
   W = mxCreateStructMatrix(1, 1, 5, fieldnames);
   plhs[0] = W;
   
   
   mxF       = mxCreateDoubleMatrix(N, N, mxREAL); F = mxGetData(mxF);
   mxdFds    = mxCreateDoubleMatrix(N, N, mxREAL); dFds = mxGetData(mxdFds);
   mxdFdsig2 = mxCreateDoubleMatrix(N, N, mxREAL); dFdsig2 = mxGetData(mxdFdsig2);
   mxdFdB    = mxCreateDoubleMatrix(N, N, mxREAL); dFdB = mxGetData(mxdFdB);
   mxdFdl    = mxCreateDoubleMatrix(N, N, mxREAL); dFdl = mxGetData(mxdFdl);
   mxdFdh    = mxCreateDoubleMatrix(N, N, mxREAL); dFdh = mxGetData(mxdFdh);
   
   F[0] = 1; F[N*N-1] = 1;
   
   /* define the slices of a gaussian with sigma2 variance */
   ndeltas = ceil(10*sqrt(sigma2)/dx); 
   if (ndeltas<50) {ndeltas = 50;}
   {
   double deltas[2*ndeltas+1];
   double ps[2*ndeltas+1];
   if (sigma2 == 0) {
       ndeltas = 1;
       deltas[0] = 0;
       ps[0] = 1;
   }
   else {
       for (i=-ndeltas; i<=ndeltas; i++){ 
           deltas[i+ndeltas]  = i * (5*sqrt(sigma2))/ndeltas;
       } 
       ndeltas = ndeltas*2+1;
       double ps_sum=0;
       for (i=0; i<ndeltas; i++){
           ps[i] = exp(-pow(deltas[i],2)/(2*sigma2));
           ps_sum += ps[i];
       }
       for (i=0; i<ndeltas; i++){ ps[i] = ps[i]/ps_sum; } 
   }
   
   /* construct F and its derivatives, stepping through the bins j where we came from: */
   double mu;
   int hp, lp;
   for (j=1; j<N-1; j++) {
       /* mu is the center of mass where the probability starting at bin j will go */
        if (fabs(lambda) < 1.0e-10){ mu = x[j] + h*dt; }
        else { mu = exp(lambda*dt)*(x[j] + h/lambda) - h/lambda; }
        
       /* now we're going to look over all the slices of the gaussian */
        for (k=0; k<ndeltas; k++){
            double s = mu + deltas[k];
            if (s <= x[0]){ 
                F[j*N] = F[j*N] + ps[k];
            }
            else if (s >= x[N-1]){ 
                F[j*N+N-1] = F[j*N+N-1] + ps[k];
            }
            else{
                /* find the bin ids whose positions are just below or just above ss[k] */
                if (x[0]<s && s<x[1]) {lp = 0; hp = 1;}
                else if (x[N-2]<s && s<x[N-1]){lp = N-2; hp = N-1;}
                else {
                    hp = ceil( (s-x[1])/dx) + 1;
                    lp = floor((s-x[1])/dx) + 1;
                }
                
                if (hp == lp) { 
                    F[j*N+lp] = F[j*N+lp] + ps[k];
                    
                    dFds[j*N+lp]   = dFds[j*N+lp]   - ps[k]/(x[lp+1]-x[lp]);
                    dFds[j*N+lp+1] = dFds[j*N+lp+1] + ps[k]/(x[lp+1]-x[lp]);
                                        
                    dFdsig2[j*N+lp]   = dFdsig2[j*N+lp]   - (ps[k]/(x[lp+1]-x[lp]))*deltas[k]/(2*sigma2);
                    dFdsig2[j*N+lp+1] = dFdsig2[j*N+lp+1] + (ps[k]/(x[lp+1]-x[lp]))*deltas[k]/(2*sigma2);
                }
                else { 
                    double dd = x[hp] - x[lp];
                    F[j*N+hp] = F[j*N+hp] + ps[k]*(s-x[lp])/dd;
                    F[j*N+lp] = F[j*N+lp] + ps[k]*(x[hp]-s)/dd;
                    
                    dFds[j*N+hp] = dFds[j*N+hp] + ps[k]/dd;
                    dFds[j*N+lp] = dFds[j*N+lp] - ps[k]/dd;
                    
                    dFdsig2[j*N+hp] = dFdsig2[j*N+hp] + (ps[k]/dd)*deltas[k]/(2*sigma2);
                    dFdsig2[j*N+lp] = dFdsig2[j*N+lp] - (ps[k]/dd)*deltas[k]/(2*sigma2);
                    
                    if (lp==0){
                        dFdB[j*N+lp] = dFdB[j*N+lp] + dB[0]*ps[k]*(x[hp]-s)/pow(dd,2);
                        dFdB[j*N+hp] = dFdB[j*N+hp] + dB[0]*ps[k]*((s-x[lp]) - dd)/pow(dd,2);        
                    }
                    if (hp==N-1){
                        dFdB[j*N+lp] = dFdB[j*N+lp] + dB[N-1]*ps[k]*(dd - (x[hp]-s))/pow(dd,2);
                        dFdB[j*N+hp] = dFdB[j*N+hp] - dB[N-1]*ps[k]*(s - x[lp])/pow(dd,2);
                    }
                }
                
            }
        }
   }
   
   
   if (fabs(lambda)<1.0e-10){
        for (i=0; i<N*N; i++){ dFdl[i] = 0;}
        for (i=0; i<N*N; i++){ dFdh[i] = dFds[i]*dt;}
   }
   else {
        double gamma = exp(lambda*dt);
        double phi = h/lambda;
        for (i=0; i<N*N; i++){
            int thisx = floor(i/N);
            dFdl[i] = dFds[i]*(x[thisx] + phi)*dt*gamma - dFds[i]*h/pow(lambda,2)*(gamma-1);
        }
        for (i=0; i<N*N; i++){ 
            dFdh[i] = dFds[i]*(gamma-1)/lambda;
        }
   }
   
   mxSetPr(mxF, F);
   
   mxSetField(W, 0, fieldnames[0], mxF);
   mxSetField(W, 0, fieldnames[1], mxdFdsig2);
   mxSetField(W, 0, fieldnames[2], mxdFdB);
   mxSetField(W, 0, fieldnames[3], mxdFdl);
   mxSetField(W, 0, fieldnames[4], mxdFdh);
   }              
}

