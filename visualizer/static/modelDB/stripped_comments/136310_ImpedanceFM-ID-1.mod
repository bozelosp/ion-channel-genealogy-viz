NEURON	{ 
	POINT_PROCESS ImpedanceFM
        POINTER vec1
        POINTER vec2
}


ASSIGNED {
	vec1
	vec2
}


VERBATIM
#include <stdlib.h>



void calc_impedances(double Zr[], double Zi[], double fmax, double df, double rext, double rmax, double dr, double R, double sigma1, double sigma2, double lambda,double epsilon, double sigmaR)
	
   {
        double epsR,sigR,w,w2,sig,eps,den,ReZ,ImZ,r,f;
	float *sigmatab;
	double PI = 3.1415927;
	int j,k,siz;


	
	
	

	
	
	siz = fmax/df + 1;
	sigmatab = (float *) malloc(sizeof(float) * siz);
	k = 0;
	for(r=rext; r<=rmax; r=r+dr) {
	  sigmatab[k] = sigma2 + (sigma1-sigma2) * exp(-(r-R)/lambda);
	  k++;
	}

	
	sigR = sigma1;
	epsR = epsilon;
	j=0;
	for(f=0; f<=fmax; f=f+df) {	
	  w = 2*PI*f;
	  w2 = w*w;
	  ReZ=0;
	  ImZ=0;
	  k=0;
	  for(r=rext; r<=rmax; r=r+dr) {	
	    
	    sig = sigmatab[k];	
	    eps = epsilon;
	    den = r*r * (sig*sig + w2 * eps*eps);
	    ReZ = ReZ + (sig*sigR + w2*eps*epsR) / den;
	    ImZ = ImZ + (sig*epsR - sigR*eps) / den;
	    k++;
	  }
	  Zr[j] = dr/(4*PI*sigmaR) * ReZ;	
	  Zi[j] = w * dr/(4*PI*sigmaR) * ImZ;
          j++;
	  
	}
	free(sigmatab);
   }


ENDVERBATIM










PROCEDURE impedance(fmax,df,rext,rmax,dr,R,sigma1,sigma2,lambda,epsilon,sigmaR) {

VERBATIM

  

  calc_impedances(&vec1,&vec2,_lfmax,_ldf,_lrext,_lrmax,_ldr,_lR,_lsigma1,\
	_lsigma2,_llambda,_lepsilon,_lsigmaR);

ENDVERBATIM
}