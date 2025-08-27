INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON { 

    
    SUFFIX EPlasInp
    RANGE C, lastrelease, lastspike, releaseat, Delay
	GLOBAL Cdur, Deadtime			 
	GLOBAL Alpha_1, Beta_1                  
	RANGE  ampa, R0_1, R1_1, Rinf_1, Rtau_1 
	GLOBAL Alpha_2, Beta_2 
    GLOBAL DurNMDA, tauNMDA
    
    RANGE  R, U, u, trec, tfac, RG 

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
 	Delay = 1.4     (ms) 

	Alpha_1 = 1.5	(/ms mM)	
	Beta_1 = 0.75	(/ms)		
	Alpha_2 = 0.25 (/ms mM) 	
	Beta_2 = 0.025 (/ms)		
    DurNMDA = 0.4               
    tauNMDA	=50                 
	tfac = 0.01     (ms)		
	trec = 500		 (ms) 	    
	U = 0.5 						
}


ASSIGNED { 

	dt		 (ms)
	v		 (mV)		
	i 		 (nA)		
	C		 (mM)		
 
	ampa
	R0_1				
	R1_1				
	Rinf_1			
	Rtau_1	(ms)		
	RG				
		 
	lastrelease	(ms)		
	lastspike	(ms) 
	releaseat

	R				
	u				
    
	
}

INITIAL {
	C = 0 
	ampa = 0
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
			Rinf_1 = C*Alpha_1 / (C*Alpha_1 + Beta_1) 
			Rtau_1 = 1 / ((Alpha_1 * C) + Beta_1) 
			R0_1 = ampa
			R_2=(1-R_2)*0.5+R_2
		}	 
						 
	} else if (q < 0) {			
		
	} else if (C > 0) {			
		C = 0. 
		R1_1 = ampa 
	} 
 	if (C > 0) {				
	   ampa = Rinf_1 + (R0_1 - Rinf_1) * exptable (- (t - lastrelease) / Rtau_1) 
	} else {					
  	   ampa = R1_1 * exptable (- Beta_1 * (t - (lastrelease + Cdur))) 
	} 

    SOLVE G_protein METHOD euler

	VERBATIM
	return 0;
	ENDVERBATIM
} 

DERIVATIVE G_protein { 							

    R_2'=-(R_2/tauNMDA)

	if (R_2<0.01) {
		RG = 0.0
	} else {
		RG = 1/( 1 + exptable(-((R_2)-DurNMDA)/0.05) ) 		
	}
	nmda' = Alpha_2 * RG * (1-nmda) - Beta_2 * nmda

    

}

FUNCTION exptable(x) { 
	TABLE  FROM -10 TO 10 WITH 2000
 
	if ((x > -10) && (x < 10)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}