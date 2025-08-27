UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX KC3
	USEION ca READ cai
	USEION k WRITE ik
	RANGE gkbar,gk,zinf,ik
}


PARAMETER {
	celsius=37	(degC)
	v		(mV)
	gkbar=.08	(mho/cm2)	
	cai = .04e-3	(mM)
	ek  = -85	(mV)
	dt		(ms)
	mon = 1
	zon = 1
}


ASSIGNED {
	ik		(mA/cm2)
	minf
	mexp
	zinf
	zexp
	gk
}

STATE {	m z }		

BREAKPOINT {
	SOLVE state

	ik = gkbar*m*z*z*(v - ek)
}






PROCEDURE state() {	
	rate(v, cai)
	m = mon * (m + mexp*(minf - m))
	z = zon * (z + zexp*(zinf - z))
	VERBATIM
	return 0;
	ENDVERBATIM
}

INITIAL {
	rate(v, cai)
	m = minf
	z = zinf
}

FUNCTION alp(v (mV), ca (mM)) (1/ms) { 
	alp = 400/(ca*1000)
}

FUNCTION bet(v (mV)) (1/ms) { 
	bet = 0.11/exp((v-35)/14.9)
}

PROCEDURE rate(v (mV), ca (mM)) { 
	LOCAL a,b
	a = alp(v,ca)
	zinf = 1/(1+a)
	zexp = (1 - exp(-dt/10))
	b = bet(v)
	minf = 7.5/(7.5+b)
	mexp = (1 - exp(-dt*(7.5+b)))
}