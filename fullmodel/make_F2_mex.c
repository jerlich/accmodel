/*=================================================================
 * make_F2_mex.c
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

/* written by BWB Aug 2011 */
#include "mex.h"
#include "matrix.h"
#include <math.h>

void
        mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
   
   int i, j, k;
   int N;
   double *x, *dB, dx;
   double sigma2;
   double lambda;
   double h;
   double dt;
   
   double *F;
   double *dFds;
   double *dFdsig2;
   double *dFdB;
   double *dFdl;
   double *dFdh;
   
   int ndeltas;
   

   
   // inputs
   N  = mxGetNumberOfElements(prhs[0]);
   x  = mxGetData(prhs[0]);
   dB = mxGetData(prhs[1]);
   sigma2 = mxGetScalar(prhs[2]);
   lambda = mxGetScalar(prhs[3]);
   h      = mxGetScalar(prhs[4]);
   dt     = mxGetScalar(prhs[5]);
   dx     = mxGetScalar(prhs[6]);
   
   
   // ouputs
   plhs[0]= mxCreateDoubleMatrix(N, N, mxREAL);
   plhs[1]= mxCreateDoubleMatrix(N, N, mxREAL);
   plhs[2]= mxCreateDoubleMatrix(N, N, mxREAL);
   plhs[3]= mxCreateDoubleMatrix(N, N, mxREAL);
   plhs[4]= mxCreateDoubleMatrix(N, N, mxREAL);
   plhs[5]= mxCreateDoubleMatrix(N, N, mxREAL);
   
   F       = mxGetData(plhs[0]);
   dFds    = mxGetData(plhs[1]);
   dFdsig2 = mxGetData(plhs[2]);
   dFdB    = mxGetData(plhs[3]);
   dFdl    = mxGetData(plhs[4]);
   dFdh    = mxGetData(plhs[5]);
   
   F[0] = 1; F[N*N-1] = 1;
   
   // define the slices of a gaussian with sigma2 variance
   if (sigma2 == 0){ ndeltas = 1; }
   else { ndeltas = ceil(10*sqrt(sigma2)/dx); }
   
   double deltas[2*ndeltas+1];
   double ps[2*ndeltas+1];
   for (i=-ndeltas; i<=ndeltas; i++){ 
       deltas[i+ndeltas]  = i * (5*sqrt(sigma2))/ndeltas;
   } 
   ndeltas = ndeltas*2+1;
   double ps_sum=0;
   for (i=0; i<ndeltas; i++){
       ps[i] = exp(-pow(deltas[i],2)/(2*sigma2));
       ps_sum += ps[i];
   }
   for (i=0; i<ndeltas; i++){ ps[i] = ps[i]/ps_sum; } // normalize
   
   // construct F and its derivatives, stepping through the bins j where we came from:
   double mu;
   int hp, lp;
   for (j=1; j<N-1; j++) {
       // mu is the center of mass where the probability starting at bin j will go
        if (abs(lambda) == 0){ mu = x[j] + h*dt; }
        else { mu = exp(lambda*dt)*(x[j] + h/lambda) - h/lambda; }
        
       // now we're going to look over all the slices of the gaussian
        for (k=0; k<ndeltas; k++){
            double s = mu + deltas[k];
            if (s <= x[0]){ // everything goes in first bin
                F[j*N] = F[j*N] + ps[k];
            }
            else if (s >= x[N-1]){ // everything goes in last bin
                F[j*N+N-1] = F[j*N+N-1] + ps[k];
            }
            else{
                // find the bin ids whose positions are just below or just
                // above ss[k]
                if (x[0]<s && s<x[1]) {lp = 0; hp = 1;}
                else if (x[N-2]<s && s<x[N-1]){lp = N-2; hp = N-1;}
                else {
                    hp = ceil( (s-x[1])/dx) + 1;
                    lp = floor((s-x[1])/dx) + 1;
                }
                
                if (hp == lp) { // if ss[k] landed exactly on a bin
                    F[j*N+lp] = F[j*N+lp] + ps[k];
                    
                    dFds[j*N+lp]   = dFds[j*N+lp]   - ps[k]/(x[lp+1]-x[lp]);
                    dFds[j*N+lp+1] = dFds[j*N+lp+1] + ps[k]/(x[lp+1]-x[lp]);
                    
                    dFdsig2[j*N+lp]   = dFdsig2[j*N+lp]   - ps[k]/(x[lp+1]-x[lp])*deltas[k]/(2*sigma2);
                    dFdsig2[j*N+lp+1] = dFdsig2[j*N+lp+1] + ps[k]/(x[lp+1]-x[lp])*deltas[k]/(2*sigma2);
                }
                else { // if ss[k] is between bins
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
   
   
   if (abs(lambda)==0){
        for (i=0; i<N*N-1; i++){ dFdl[i] = 0;}
        for (i=0; i<N*N-1; i++){ dFdh[i] = dFds[i]*dt;}
   }
   else {
        double gamma = exp(lambda*dt);
        double phi = h/lambda;
        for (i=0; i<N*N-1; i++){
            int thisx = floor(i/N);
            dFdl[i] = dFds[i]*(x[thisx] + phi)*dt*gamma - dFds[i]*h/pow(lambda,2)*(gamma-1);
        }
        for (i=0; i<N*N-1; i++){ 
            dFdh[i] = dFds[i]*(gamma-1)/lambda;
        }
   }
   
                  
}

