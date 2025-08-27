INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAaSynapse
	USEION cl READ ecl VALENCE 1
	
	RANGE onset, gmaxIPSP, e, g, i, w
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) 	= (nanoamp)
	(mV)	= (millivolt)
	(nS) 	= (nanomho)
}

PARAMETER {
	onset= 25	(ms)
	gmaxIPSP= 0	(nS)
	w= 1 				
	ecl		(mV)
	v		(mV)
	celsius		(degC)
}

ASSIGNED { 
	i 		(nA)  
	g 		(nS)
	tadj
}

UNITSOFF
INITIAL {
        tadj = 3^((celsius-23.5)/10)
}    

BREAKPOINT { LOCAL tt
	tt= (t-onset)*tadj
	if ((t>onset)&&(tt<740)) {
	
	  g = w*gmaxIPSP * exp(-tt/25) * (1-exp(-tt/1.0))/0.84
	}
	else {g = 0}
	
	i = g * (v-(-ecl))	
}
UNITSON