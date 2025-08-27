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
	SUFFIX mykca
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gkbar, ik
	GLOBAL oinf, tau
}

UNITS {
	FARADAY = (faraday)  (kilocoulombs)
	R = 8.313424 (joule/degC)
}

PARAMETER {
	v		(mV)
	dt		(ms)
	ek		(mV)
	celsius = 20	(degC)
	gkbar = 0.01	(mho/cm2)	
	cai = 1e-3	(mM)
	d1 = 0.84
	d2 = 1.0
	k1 = 0.18	(mM)
	k2 = 0.011	(mM)
	bbar = 0.28	(/ms)
	abar = 0.48	(/ms)
}


ASSIGNED {
	ik		(mA/cm2)
	oinf
	tau		(ms)
}

STATE {	o }		

BREAKPOINT {
	SOLVE state
	ik = gkbar*o*(v - ek) 
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
	alp = abar/(1 + exp1(k1,d1,v)/ca)
}

FUNCTION bet(v (mV), ca (mM)) (1/ms) { 
	bet = bbar/(1 + ca/exp1(k2,d2,v))
}  

FUNCTION exp1(k (mM), d, v (mV)) (mM) { 
	exp1 = k*exp(-2*d*FARADAY*v/R/(273.15 + celsius))
}

PROCEDURE rate(v (mV), ca (mM)) { 
	LOCAL a
	a = alp(v,ca)
	tau = 1/(a + bet(v, ca)) 
	oinf = a*tau             
	fac = (1 - exp(-dt/tau))
}