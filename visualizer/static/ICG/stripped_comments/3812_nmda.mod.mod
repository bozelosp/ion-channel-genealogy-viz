INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS NMDA_muki
	RANGE   g, gmax, e, i, tau1, tau2, gama, ni, Mgc, onset 
	NONSPECIFIC_CURRENT i
	 
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	onset = 0 (ms)		
  	gmax = 0.001 (umho) 	
	tau1 = 40 (ms)
	tau2 = 0.33 (ms)
	e = 0 (mV)
	gama = 0.06 (/mv)
	ni = 0.33 (/mM)
	Mgc = 1 (mM)	
}


ASSIGNED {
	v  (mV)		
	i  (nA)		
	g  (umho)		
}


 
BREAKPOINT {
	LOCAL  td
	if ( t < onset ) { 
		i = 0
	}  else {
	        td = t - onset
		g = calc_nmda ( td ) 
		i = g*(v - e)
	}
}







FUNCTION calc_nmda ( td ) {
	LOCAL gama_exp

	if ( gama * v < -14 ) {  		 
		gama_exp = 1e6
	} else {
		gama_exp = exp(-gama*v)	  	 
	}
	calc_nmda =  gmax * ( exp(-td/tau1) - exp(-td/tau2) )/ (1 + ni * Mgc * gama_exp )
	
}