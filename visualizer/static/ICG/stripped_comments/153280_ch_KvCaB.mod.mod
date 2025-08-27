VERBATIM
#include <stdlib.h> 
ENDVERBATIM

UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}

NEURON {
	SUFFIX ch_KvCaB
	USEION k READ ek WRITE ik
	USEION ca READ cai 
	RANGE gmax, g, ik
	RANGE myi
	GLOBAL oinf, otau	
    THREADSAFE
}

UNITS {
	FARADAY = (faraday)  (kilocoulombs)
	R = 8.313424 (joule/degC)
}

PARAMETER {	
	gmax=.01	(mho/cm2)	

	d1 = .84
	d2 = 1.	
	k1 = .48e-3	(mM)
	k2 = .13e-6	(mM)
	cai = 5.e-5	(mM)
	myshift = 25 (mV)
	
	abar = .28	(/ms)
	bbar = .48	(/ms)
	
	st=1		(1)
}

ASSIGNED {	
      celsius (degC) 
	v			(mV)


	ek			(mV)
	ik			(mA/cm2)

	oinf
	otau		(ms)
	g		(mho/cm2)
	myi (mA/cm2)
}

INITIAL {
        rate(v,cai)
        o=oinf
}

STATE {	o }		

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gmax*o^st
	ik = g*(v - ek)
	myi = ik
}

DERIVATIVE state {	
	rate(v, cai)
	o' = (oinf - o)/otau
}

FUNCTION alp(v (mV), c (mM)) (1/ms) { 
	alp = c*abar/(c + exp1(k1,d1,v))
}

FUNCTION bet(v (mV), c (mM)) (1/ms) { 
	bet = bbar/(1 + c/exp1(k2,d2,v))
}

FUNCTION exp1(k (mM), d, v (mV)) (mM) { 
	exp1 = k*exp(-2*d*FARADAY*(v-myshift)/R/(273.15 + celsius))
}

PROCEDURE rate(v (mV), c (mM)) { 
	LOCAL a
	a = alp(v,c)
	otau = 1/(a + bet(v, c))
	oinf = a*otau
}