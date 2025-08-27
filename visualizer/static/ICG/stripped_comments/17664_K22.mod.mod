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
	SUFFIX K22
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gkbar,gk,zinf,ik
}


PARAMETER {
	celsius=37	(degC)
	v		(mV)
	gkbar=.00039	(mho/cm2)	
	cai = .04e-3	(mM)
	
	dt		(ms)
}


ASSIGNED {
        ek (mV)
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
	m = m + mexp*(minf - m)
	z = z + zexp*(zinf - z)
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
	alp = 2/(ca*1000)
}

FUNCTION bet(v (mV)) (1/ms) { 
	bet = 0.075/exp((v+5)/10)
}

PROCEDURE rate(v (mV), ca (mM)) { 
	LOCAL a,b
	a = alp(v,ca)
	zinf = 1/(1+a)
	zexp = (1 - exp(-dt/10))
	b = bet(v)
	minf = 25/(25+b)
	mexp = (1 - exp(-dt*(25+b)))
}