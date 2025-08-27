INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON { 

   SUFFIX IPlasSom
	RANGE C, lastrelease, lastspike, releaseat, Delay
	GLOBAL Cdur, Deadtime		
	GLOBAL Alpha, Beta 
	RANGE  gaba, R0, R1, Rinf, Rtau 
   
   GLOBAL U, trec, tfac
   RANGE  R, u, RG 

} 
 
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

STATE {
	nmda				
	R_2
}

PARAMETER {

	Cdur	= 1	  (ms)		
	Deadtime = 1  (ms)		
	Delay = 0.6

	Alpha = 0.5		(/ms mM)	
	Beta = 0.25		(/ms)		
 
	
	tfac = 20     (ms)		
	trec = 700		 (ms) 	
	U = 0.25 						
	
}


ASSIGNED { 
 
	dt		 (ms)
	v		 (mV)		
	i 		 (nA)		
	C		 (mM)		
 
	gaba
	R0				
	R1				
	Rinf			
	Rtau	(ms)		
		 
	lastrelease	(ms)		
	lastspike	(ms) 
	releaseat

	R				
	u				
	
}

INITIAL {
	C = 0
	R1 = 0

	gaba = 0
	lastrelease = -9e4 
	lastspike   = -9e4 
	releaseat   = -9e4 
 
	R = 1 
	u = U 
 
}

BREAKPOINT {
    SOLVE release

}

PROCEDURE release() { LOCAL q
    
	
	q = (t - lastspike)		
						 
	if (q > Deadtime) {		
		if (v > 0) {		
			lastspike = t 
			releaseat = t + Delay 
		} 
	} 
 	
	q = (t - lastrelease -Cdur)			
	if (q > Deadtime) {          
		if (t > releaseat - dt/2 && t < releaseat + dt/2) { 
			lastrelease = t 
     		u = u*(exptable(-q/tfac)) + U*(1-u*exptable(-q/tfac)) 
			R = R*(1 - u)*exptable(-q/trec) + 1 - exptable(-q/trec)		 
			C = R*u				
			Rinf = C*Alpha / (C*Alpha + Beta) 
			Rtau = 1 / ((Alpha * C) + Beta) 
			R0 = gaba
		}	 						 
	} else if (q < 0) {			
		
	} else if (C > 0) {			
		C = 0. 
		R1 = gaba 
	} 

	if (C > 0) {				
	   gaba = Rinf + (R0 - Rinf) * exptable (- (t - lastrelease) / Rtau) 
	} else {					
  	   gaba = R1 * exptable (- Beta * (t - (lastrelease + Cdur))) 
	} 

	VERBATIM
	return 0;
	ENDVERBATIM
} 


FUNCTION exptable(x) { 
	TABLE  FROM -10 TO 10 WITH 2000
 
	if ((x > -10) && (x < 10)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}