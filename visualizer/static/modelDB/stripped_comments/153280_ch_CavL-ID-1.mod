VERBATIM
#include <stdlib.h> 
ENDVERBATIM

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
}

NEURON {
	SUFFIX ch_CavL
	USEION ca READ cai, cao, eca WRITE ica VALENCE 2 
    RANGE gmax, g, cai, ica, eca
 	RANGE myi
    RANGE minf, mtau
    THREADSAFE
}

PARAMETER {
	v (mV)
    celsius (degC) 
	gmax		 (mho/cm2)
	ki=.001 (mM)
	cai (mM)
	cao (mM)
	tfa=1
}

STATE {
	m
}

ASSIGNED {
	g (mho/cm2)
	minf
	mtau   (ms)
	myi (mA/cm2)
	ica (mA/cm2)
	eca (mV)   
}

INITIAL {
	rate(v)
	m = minf
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = gmax*m*m*h2(cai)
	ica = g*ghk(v,cai,cao)
	myi = ica
}

FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
	LOCAL f
	f = KTF(celsius)/2
	ghk=-f*(1. - (ci/co)*exp(v/f))*vtrap(v,f)/f
}

FUNCTION KTF(celsius (DegC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION vtrap(x,y) {  
	if (fabs(x/y) < 1e-6) {
			vtrap = y*(1 - x/y/2)
	}else{  
			vtrap = x/(exp(x/y) - 1)
	}
}

FUNCTION alp(v(mV)) (1/ms) {
	alp = 15.69*vtrap((-1.0*v+81.5),10.0)
}

FUNCTION bet(v(mV)) (1/ms) {
	bet = 0.29*exp(-v/10.86)
}

DERIVATIVE state {  
        rate(v)
        m' = (minf - m)/mtau
}

PROCEDURE rate(v (mV)) { 
        LOCAL a
        a = alp(v)
        mtau = 1/(tfa*(a + bet(v)))
        minf = tfa*a*mtau
}