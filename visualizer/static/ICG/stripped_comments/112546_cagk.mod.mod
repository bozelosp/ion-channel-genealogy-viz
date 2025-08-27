UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


NEURON {
	SUFFIX cagk
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gbar,gkca,ik
	GLOBAL oinf, tau
}

UNITS {
	FARADAY = (faraday)  (kilocoulombs)
	R = 8.313424 (joule/degC)
}

PARAMETER {
	celsius		(degC)
	v		(mV)
	gbar=.01	(mho/cm2)	
	cai=50e-6 		(mM)	
	ek		(mV)

	d1 = .84
	d2 = 1.
	k1 = .48e-3	(mM)
	k2 = .13e-6	(mM)
	abar = .28	(/ms)
	bbar = .48	(/ms)
        st=1            (1)
}

ASSIGNED {
	ik		(mA/cm2)
	oinf
	tau		(ms)
        gkca          (mho/cm2)
}

INITIAL {
        rate(v,cai)
        o=oinf
}

STATE {	o }		

BREAKPOINT {
	SOLVE state METHOD cnexp
	gkca = gbar*o^st
	ik = gkca*(v - ek)
}

DERIVATIVE state {	
	rate(v, cai)
	o' = (oinf - o)/tau
}

FUNCTION alp(v (mV), c (mM)) (1/ms) { 
	alp = c*abar/(c + exp1(k1,d1,v))
}

FUNCTION bet(v (mV), c (mM)) (1/ms) { 
	bet = bbar/(1 + c/exp1(k2,d2,v))
}

FUNCTION exp1(k (mM), d, v (mV)) (mM) { 
	exp1 = k*exp(-2*d*FARADAY*v/R/(273.15 + celsius))
}

PROCEDURE rate(v (mV), c (mM)) { 
	LOCAL a
	a = alp(v,c)
	tau = 1/(a + bet(v, c))
	oinf = a*tau
	if (oinf<1.e-6) {oinf=0}
}