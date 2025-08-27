VERBATIM
#include <stdlib.h> 
ENDVERBATIM

NEURON {
SUFFIX iconc_Ca
USEION ca READ cai, ica, eca WRITE eca, cai VALENCE 2
RANGE caiinf, catau, cai, eca
THREADSAFE
}

UNITS {
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (milli/liter)
	(mA) = (milliamp)
	FARADAY = 96520 (coul)
	R = 8.3134	(joule/degC)
}

INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

PARAMETER {
    celsius (degC) 
	depth = 200 (nm)	
	catau = 9 (ms)
	caiinf = 50.e-6 (mM)	
			
			
	cao = 2 (mM)
	ica (mA/cm2)
}

ASSIGNED {
	eca (mV)
}

STATE {
	cai
}


INITIAL {
	
	cai = caiinf	
	eca = ktf() * log(cao/caiinf)	
}


BREAKPOINT {
	SOLVE integrate METHOD derivimplicit
	eca = ktf() * log(cao/cai)	
}

DERIVATIVE integrate {
cai' = -(ica)/depth/FARADAY * (1e7) + (caiinf - cai)/catau
}

FUNCTION ktf() (mV) {
	ktf = (1000)*R*(celsius +273.15)/(2*FARADAY)
}