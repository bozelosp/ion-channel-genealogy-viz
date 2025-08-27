VERBATIM
#include <stdlib.h> 
ENDVERBATIM

UNITS {
        (molar) = (1/liter)
        (mM)    = (millimolar)
	(mA)	= (milliamp)
	(mV)	= (millivolt)
}

NEURON {
	SUFFIX ch_KCaS
	USEION k READ ek WRITE ik 
	USEION ca READ cai 
	RANGE g, gmax, qinf, qtau, ik
	RANGE myi
    THREADSAFE
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
      celsius (degC) 
	v		(mV)
	dt		(ms)
	gmax = 1.0  (mho/cm2)
	ek	(mV)
	cai (mM)
}

STATE { q }

ASSIGNED {
	ik (mA/cm2) g(mho/cm2) qinf qtau (ms) qexp
	myi (mA/cm2)
}


BREAKPOINT {          
	SOLVE state
        g = gmax * q*q
	ik = g * (v-ek)
	myi = ik
}

UNITSOFF

INITIAL {
	rate(cai)
	q=qinf
}

PROCEDURE state() {  
	q = q + (qinf-q) * qexp
	rate(cai)
}

LOCAL q10
PROCEDURE rate(cai) {  
	LOCAL alpha, beta, tinc
	
	q10 = 3^((celsius - 34)/10) 
		
alpha = 1.25e1 * cai * cai
beta = 0.00025 





	qtau = 1 /(alpha + beta)/q10
	qinf = alpha * qtau
	tinc = -dt
	qexp = 1 - exp(tinc/qtau)
}

UNITSON