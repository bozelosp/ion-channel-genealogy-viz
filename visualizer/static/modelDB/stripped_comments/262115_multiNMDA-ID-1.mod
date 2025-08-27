DEFINE MAXSYNNMDA 250
VERBATIM
static int MAXSYNNMDA = 250;
ENDVERBATIM

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS multiNMDA
	NONSPECIFIC_CURRENT i
	RANGE Ron, Roff, ri, nsyn, non, g, gmax, B
	GLOBAL Cmax, Cdur, Alpha, Beta, Erev
	GLOBAL Prethresh, Deadtime, Rinf, Rtau, mg
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt		(ms)

	Cmax	= 0.5	(mM)		
	Cdur	= 0.3	(ms)		
	Alpha	= 0.11	(/ms mM)	
	Beta	= 0.0064 (/ms)		
	Erev	= 0	(mV)		
	gmax		(umho)		

	Prethresh = 0 			
	Deadtime = 1	(ms)		
	mg	= 1    (mM)		
}


ASSIGNED {
	on[MAXSYNNMDA]			
	TL[MAXSYNNMDA]	(ms)		
	ri[MAXSYNNMDA]			
	lastrelease[MAXSYNNMDA] (ms)	
	Ron				
	Roff				
	nsyn				
	non				

	Rinf				
	Rtau		(ms)		

	v		(mV)		
	i 		(nA)		
	g 		(umho)		

	trel		(ms)		
	ptr_array_nmda			
	B				
}

INITIAL { LOCAL j
	FROM j=0 TO nsyn-1 {
		on[j] = 0
		TL[j] = -9e9
		lastrelease[j] = -9e9
		ri[j] = 0
	}
	Ron = 0
	Roff = 0
	non = 0

	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / ((Alpha * Cmax) + Beta)
}

BREAKPOINT {
   if(gmax > 0) {
	SOLVE release
	B = mgblock(v)		
	g = gmax * (Ron+Roff) * B
	i = g*(v - Erev)
   } else {
	i = 0
   }
}

PROCEDURE release() { LOCAL q,j

  FROM j=0 TO nsyn-1 {	

    trel = ((t - lastrelease[j]) - Cdur)	

    if (trel > Deadtime) {			
				
	if (presynaptic(j) > Prethresh) {	
	  on[j] = 1			
	  non = non + 1
	  lastrelease[j] = t		
	  ri[j] = ri[j] * exptable( - Beta * (t-TL[j]))
					
	  TL[j] = t			
	  Ron = Ron + ri[j]		
	  Roff = Roff - ri[j]		
	  if(Roff < 1e-9) { Roff = 0 }	
	}
						
    } else if (trel < 0) {			

		
	
    } else if (on[j] > 0) {			
	on[j] = 0			
	non = non - 1
	ri[j] = Rinf + (ri[j]-Rinf) * exptable(- (t-TL[j]) / Rtau)
					
	TL[j] = t			
	Ron = Ron - ri[j]		
	Roff = Roff + ri[j]		
	if(Ron < 1e-9) { Ron = 0 }		
    }

  }


  if(Roff > 0) {			
     Roff = Roff * exptable(- Beta * dt)
     if(Roff < 1e-9) { Roff = 0 }	
  }

  if(non > 0) {				
    q = non * Rinf
    Ron = q + (Ron - q) * exptable(- dt / Rtau) 
  }

}


FUNCTION exptable(x) { 
	TABLE  FROM -25 TO 25 WITH 10000

	if ((x > -25) && (x < 25)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}

























FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}











VERBATIM
#define ppnmda ((double***)(&(ptr_array_nmda)))
extern double* hoc_pgetarg();
ENDVERBATIM





PROCEDURE allocate(n) {
  VERBATIM
	if (*ppnmda) {
	   free(*ppnmda);
	}
	*ppnmda = ((double**) hoc_Ecalloc((int)_ln, sizeof(double *))), hoc_malchk();
  ENDVERBATIM
}





FUNCTION presynaptic(index) {
  VERBATIM
	if(_lindex >= nsyn) {
	   printf("Warning: attempt to use pointer outside range\n");
	   printf(" trying to use pointer number %d\n",(int)_lindex);
	   printf(" but number of defined pointers was nsyn=%d.\n",(int) nsyn);
	}
	_lpresynaptic = *((*ppnmda)[(int)_lindex]);
  ENDVERBATIM
}







PROCEDURE addlink() {
  VERBATIM
	if(++nsyn > MAXSYNNMDA) {
	  printf("Exceeding maximum of allowed links MAXSYNNMDA=%d\n",MAXSYNNMDA);
	  printf("  edit the nmodl code to increase the maximum allowed.\n");
	  exit(-1);
	}
	(*ppnmda)[(int)(nsyn-1)] = hoc_pgetarg(1);
  ENDVERBATIM
}