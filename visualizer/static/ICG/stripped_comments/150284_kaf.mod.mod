UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

NEURON {
	SUFFIX kaf
	USEION k READ ek WRITE ik
	RANGE gkbar, ik, mshift, hshift
}

PARAMETER {
	gkbar = 0.21	(S/cm2)		

	mvhalf = -10.0	(mV)		
	mslope = -17.7	(mV)		
	mshift = 0	(mV)
	
	hvhalf = -75.6	(mV)		
	hslope	= 10	(mV)		
	hshift = 0	(mV)
	htau = 14	(ms)			

	qfact = 3
	power = 2
}

ASSIGNED { 
	v 		(mV)
    ik 		(mA/cm2)
	ek 		(mV)

	minf
	hinf
}

STATE {
    m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	ik  = gkbar * m^power * h * (v-ek)
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

FUNCTION_TABLE mtau (v(mV))  (ms)	

DERIVATIVE states {  
	rates(v)
	m' = (minf - m) / (mtau(v) / qfact)
	h' = (hinf - h) / (htau / qfact)
}


PROCEDURE rates( v(mV) ) {  
	TABLE minf, hinf DEPEND mshift, hshift, hslope
		FROM -200 TO 200 WITH 201
			minf = 1  /  ( 1 + exp( (v - mvhalf - mshift) / mslope) ) 
			hinf = 1  /  ( 1 + exp( (v - hvhalf - hshift) / hslope) ) 
}