UNITS {
	
	(mA) = (milliamp)
	(mV) = (millivolt)
	(pS) = (picosiemens)
	(um) = (micron)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {	
	
	dt	(ms)
	v 	(mV)
	celsius	(degC)
	Ca	=	0.0388
	Cb	= 	0.00168
	za	=	0.90979
	zb	=	1.23645
	gbar	= 	20 	(pS/um2)	
	temp	= 	34	(degC)		
	q10  	= 	3.0				
}

NEURON {
	
	SUFFIX kv7
	USEION k READ ek WRITE ik
	RANGE gbar, ik
	THREADSAFE
}

STATE { m }

ASSIGNED {
	
	ik (mA/cm2)
	gk (pS/um2)
	ek (mV)
	tadj
}

INITIAL {
	
	m=alpha(v)/(beta(v)+alpha(v))
	tadj = q10^((celsius - temp)/10)
}

BREAKPOINT {
	
	SOLVE state METHOD cnexp
	tadj = q10^((celsius - temp)/10)	
	ik = (1e-4) * gbar * m * (v-ek)
}

FUNCTION alpha(v(mV)) {
	
	alpha = tadj*Ca*exp(za*v*0.037788)
}

FUNCTION beta(v(mV)) {
	
	beta = tadj*Cb*exp(-zb*v*0.037788)
}

DERIVATIVE state {
	
	m' = (1-m)*alpha(v) - m*beta(v)
}