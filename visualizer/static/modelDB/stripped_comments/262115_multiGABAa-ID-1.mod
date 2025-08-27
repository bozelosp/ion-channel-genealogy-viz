DEFINE MAXSYNGABAA 260
VERBATIM
static int MAXSYNGABAA = 260;
ENDVERBATIM

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS multiGABAa
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
	Alpha	= 20	(/ms mM)	
	Beta	= 0.162	(/ms)		
	Erev	= -80	(mV)		
	gmax		(umho)		

	Prethresh = 0 			
	Deadtime = 1	(ms)		
}


ASSIGNED {
	on[MAXSYNGABAA]			
	TL[MAXSYNGABAA]	(ms)		
	ri[MAXSYNGABAA]			
	lastrelease[MAXSYNGABAA] (ms)	
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
	ptr_array_gabaa			
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
#define ppgabaa ((double***)(&(ptr_array_gabaa)))
extern double* hoc_pgetarg();
ENDVERBATIM





PROCEDURE allocate(n) {
  VERBATIM
	if (*ppgabaa) {
	   free(*ppgabaa);
	}
	*ppgabaa = ((double**) hoc_Ecalloc((int)_ln, sizeof(double *))), hoc_malchk();
  ENDVERBATIM
}





FUNCTION presynaptic(index) {
  VERBATIM
	if(_lindex >= nsyn) {
	   printf("Warning: attempt to use pointer outside range\n");
	   printf(" trying to use pointer number %d\n",(int)_lindex);
	   printf(" but number of defined pointers was nsyn=%d.\n",(int) nsyn);
	}
	_lpresynaptic = *((*ppgabaa)[(int)_lindex]);
  ENDVERBATIM
}







PROCEDURE addlink() {
  VERBATIM
	if(++nsyn > MAXSYNGABAA) {
	  printf("Exceeding maximum of allowed links MAXSYNGABAA=%d\n",MAXSYNGABAA);
	  printf("  edit the nmodl code to increase the maximum allowed.\n");
	  exit(-1);
	}
	(*ppgabaa)[(int)(nsyn-1)] = hoc_pgetarg(1);
  ENDVERBATIM
}