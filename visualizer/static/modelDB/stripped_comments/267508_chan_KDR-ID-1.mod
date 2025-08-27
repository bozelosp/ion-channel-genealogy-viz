INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

UNITS { 
	(mV) = (millivolt) 
	(mA) = (milliamp) 
} 
NEURON { 
	SUFFIX KDR
	USEION k READ ek WRITE ik
	GLOBAL gmax, minf, mtau
	RANGE ik
}
PARAMETER { 
	gmax = 0 	(mho/cm2)
} 
ASSIGNED { 
	v 		(mV)
	ek 		(mV) 
	ik 		(mA/cm2) 
	minf 		(1)
	mtau 	(ms)	
} 
STATE {
	m 
}
BREAKPOINT { 
	SOLVE states METHOD cnexp
	ik = gmax * m * m * m * m * (v - ek)
} 
INITIAL { 
	settables(v) 
	m = minf	
} 
DERIVATIVE states { 
	settables(v)
		m' = ( minf - m ) / mtau 
}		


PROCEDURE settables(v(mV)) {LOCAL alpham, betam
	TABLE minf, mtau FROM -120 TO 40 WITH 161
	alpham	= -(0.616 + (0.014 * v))/(exp(-( 44 + v ) / 2.3 ) - 1)
	betam	= 0.0043 / exp((44 + v)/ 34 )
	minf  		= alpham / (alpham + betam)
	mtau 	= 1 / (alpham + betam)
}