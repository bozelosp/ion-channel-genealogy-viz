UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 100 (ms)}

NEURON {
	SUFFIX Ic
	USEION ca READ cai
	
	USEION k READ ek WRITE ik
        RANGE gbar
	GLOBAL oinf, tau
}

UNITS {
}

PARAMETER {
	v		(mV)
	gbar=.0873	(mho/cm2)	
	
	
	dt		(ms)
	

	f = 0.0851
	g = 0.077
	k1 = 1.5e-3	(mM)
	k2 = 1.5e-4	(mM)
	bbar = 1.5	(/ms)
	abar = 2.5	(/ms)
}


ASSIGNED {
        cai (mM)
        ek (mV)
	ik 		(mA/cm2)
	oinf
	tau		(ms)
}

STATE {	o }		

BREAKPOINT {
	SOLVE state
	ik = gbar*o*(v - ek)
}

LOCAL fac




PROCEDURE state() {	
	rate(v, cai)
	o = o + fac*(oinf - o)
	VERBATIM
	return 0;
	ENDVERBATIM
}

INITIAL {
	rate(v, cai)
	o = oinf
}

FUNCTION alp(v (mV), ca (mM)) (1/ms) { 
	alp = abar/(1 + exp1(k1,f,v)/ca)
}

FUNCTION bet(v (mV), ca (mM)) (1/ms) { 
	bet = bbar/(1 + ca/exp1(k2,g,v))
}

FUNCTION exp1(k (mM), d, v (mV)) (mM) { 
	exp1 = k*exp(-d*v)
}

PROCEDURE rate(v (mV), ca (mM)) { 
	LOCAL a
	a = alp(v,ca)
	tau = 1/(a + bet(v, ca))
	oinf = a*tau
	fac = (1 - exp(-dt/tau))
}