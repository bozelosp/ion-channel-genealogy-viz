NEURON {
	SUFFIX hBKKCa
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE  gkbar,ik
}

UNITS {
	(molar) = (1/liter)
	(mM)	= (millimolar)
	(S)  	= (siemens)
	(mA) 	= (milliamp)
	(mV) 	= (millivolt)
}

PARAMETER {
	gkbar	= 0.12 (S/cm2)
	
	mvhalf = 2.6	(mV)		
	mslope = -17.0	(mV)		
	mshift = 0	(mV)
	
	hvhalf = -86.15	(mV)		
	hslope = 11.61	(mV)		
	hshift = 0	(mV)
	
	pvhalf = -7.4 	(mM)		
	pslope = -0.65	(mM)		
	pshift = 0		(mM)

	mtau = 0.133	(ms)		
	htau = 1.27	(ms)		

	mpower = 3
	ppower = 8					
	hpower = 1
}

ASSIGNED {
    v       (mV)
    cai	(mM)
	ik	(mA/cm2)
	ek 	(mV)
	po	
	minf
	hinf
}

STATE { m h }

BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = po^ppower * gkbar * m^mpower * h^hpower * (v - ek)
}

DERIVATIVE state {
	calc_po(cai)
	settables(v)
	m' = (minf - m) / mtau
	h' = (hinf - h) / htau
}

PROCEDURE settables( v (mV) ) {
	TABLE minf, hinf DEPEND mshift, hshift
        FROM -100 TO 100 WITH 201
		minf = 1  /  ( 1 + exp( (v-mvhalf-mshift) / mslope) )
		hinf = 1  /  ( 1 + exp( (v-hvhalf-hshift) / hslope) ) 
}

PROCEDURE calc_po( cai(mM) ) {
		po = 1  /  ( 1 + exp( (cai-pvhalf-pshift) / pslope) ) 
}