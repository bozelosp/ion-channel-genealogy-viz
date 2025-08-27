UNITS {
        (molar) = (1/liter)
        (mM)    = (millimolar)
	(mA)	= (milliamp)
	(mV)	= (millivolt)
}

NEURON {
	SUFFIX gskch
	
	
	
	
	USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE gsk, gskbar, qinf, qtau
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	celsius=6.3 (degC)
	v		(mV)
	dt		(ms)
	gskbar = 1.0 (mho/cm2)
	
	
	
	
	
}

STATE { q }

ASSIGNED {
        ek (mV)
        cai (mM)
	ik (mA/cm2) gsk (mho/cm2) qinf qtau (ms) qexp
}


BREAKPOINT {          
	SOLVE state
        gsk = gskbar * q*q
	ik = gsk * (v-ek)
}

UNITSOFF

INITIAL {
	
	q=qinf
	rate(cai)
	
}


PROCEDURE state() {  
	
	rate(cai)
	q = q + (qinf-q) * qexp
	VERBATIM
	return 0;
	ENDVERBATIM
}

LOCAL q10
PROCEDURE rate(cai) {  
	LOCAL alpha, beta, tinc
	q10 = 3^((celsius - 6.3)/10)
		
alpha = 1.25e1 * cai * cai
beta = 0.00025 





	qtau = 1 / (alpha + beta)
	qinf = alpha * qtau
	tinc = -dt*q10
	qexp = 1 - exp(tinc/qtau)*q10
}

UNITSON