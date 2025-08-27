UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
	FARADAY = (faraday)  (kilocoulombs)
	R = (k-mole) (joule/degC)
}


NEURON {
	SUFFIX skkca
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gkbar, ik, qfact, abar, bbar, stau
	GLOBAL oinf, tau
}

PARAMETER {
	stau = 1
	qfact = 0.72
	celsius_sk	= 22	(degC) 
	v		(mV)
	gkbar=0.175	(mho/cm2)	
	cai		(mM) 
	ek		(mV)

	d1 = .84	      
	d2 = 1.0			
	k1 = .18	(mM)
	k2 = .011	(mM)
	abar = .48	(/ms)
	bbar = .28	(/ms) 
}

ASSIGNED {
	ik		(mA/cm2)
	oinf
	tau		(ms)
}

STATE {	o }		

BREAKPOINT {
	SOLVE state METHOD cnexp
	ik = gkbar*o*(v - ek)
}

DERIVATIVE state {
	rate(v, cai)
	o' = (oinf - o)/(tau/qfact)
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
	exp1 = k*exp(-2*d*FARADAY*v/R/(273.15 + celsius_sk))
}

PROCEDURE rate(v (mV), ca (mM)) { 
	LOCAL a
	a = alp(v,ca)
	tau = stau/(a + bet(v, ca))
	oinf = a*tau
}