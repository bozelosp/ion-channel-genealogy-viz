DEFINE MAXSYNAMPA 250
VERBATIM
static int MAXSYNAMPA = 250;
ENDVERBATIM

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS multiAMPA
	NONSPECIFIC_CURRENT i
	RANGE Ron, Roff, ri, nsyn, non, g, gmax
	GLOBAL Cmax, Cdur, Alpha, Beta, Erev, Prethresh, Deadtime, Rinf, Rtau
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
	Alpha	= 0.94	(/ms mM)	
	Beta	= 0.18	(/ms)		
	Erev	= 0	(mV)		
	gmax		(umho)		

	Prethresh = 0 			
	Deadtime = 1	(ms)		
}


ASSIGNED {
	on[MAXSYNAMPA]			
	TL[MAXSYNAMPA]	(ms)		
	ri[MAXSYNAMPA]			
	lastrelease[MAXSYNAMPA] (ms)	
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
	ptr_array_ampa			
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
	g = gmax * (Ron+Roff)
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
































VERBATIM
#define ppampa ((double***)(&(ptr_array_ampa)))
extern double* hoc_pgetarg();
ENDVERBATIM





PROCEDURE allocate(n) {
  VERBATIM
	if (*ppampa) {
	   free(*ppampa);
	}
	*ppampa = ((double**) hoc_Ecalloc((int)_ln, sizeof(double *))), hoc_malchk();
  ENDVERBATIM
}





FUNCTION presynaptic(index) {
  VERBATIM
	if(_lindex >= nsyn) {
	   printf("Warning: attempt to use pointer outside range\n");
	   printf(" trying to use pointer number %d\n",(int)_lindex);
	   printf(" but number of defined pointers was nsyn=%d.\n",(int) nsyn);
	}
	_lpresynaptic = *((*ppampa)[(int)_lindex]);
  ENDVERBATIM
}







PROCEDURE addlink() {
  VERBATIM
	if(++nsyn > MAXSYNAMPA) {
	  printf("Exceeding maximum of allowed links MAXSYNAMPA=%d\n",MAXSYNAMPA);
	  printf("  edit the nmodl code to increase the maximum allowed.\n");
	  exit(-1);
	}
	(*ppampa)[(int)(nsyn-1)] = hoc_pgetarg(1);
  ENDVERBATIM
}