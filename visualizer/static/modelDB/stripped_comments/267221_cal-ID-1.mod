NEURON {
	SUFFIX cal
	USEION ca READ cai, cao WRITE ica
	RANGE gcalbar, ica
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {	
	gcalbar = 0 (S/cm2)	
	ki = 0.001 (mM)
	tfa = 5				
}

ASSIGNED { 
	v (mV)
	eca (mV)		
	cai (mM)		
	cao (mM)		
	celsius (degC)	

	ica (mA/cm2)
	gcal (S/cm2)
	minf (1)
	taum (ms)
}

STATE {	
	m
} 

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcal = gcalbar*m*h2(cai)	
	ica = gcal*ghk(v, cai, cao)	
}

DERIVATIVE states {
	rates (v)
	m' = (minf - m)/taum
}

INITIAL {
	rates(v)
	m = minf
}

FUNCTION h2(cai (mM)) {
	h2 = ki/(ki+cai)
}

FUNCTION ghk(v (mV), ci (mM), co (mM)) (mV) {
	LOCAL nu, f
	f = KTF(celsius)/2
	nu = v/f
	ghk = -f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {	
	KTF = (0.0853(mV/degC)*(celsius + 273.15(degC)))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	} else {
		efun = z/(exp(z) - 1)
	}
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {  
		vtrap = 1(/mV)*x/(exp(x/y) - 1)
	}
}

FUNCTION alpm(v (mV)) (/ms){
	alpm = 0.055(/ms)*vtrap(-(v+27.01(mV)), 3.8(mV))
}

FUNCTION betm(v (mV)) (/ms){
	betm =0.94(/ms)*exp(-(v + 63.01(mV))/17(mV))
}

PROCEDURE rates(v (mV)) { 
	taum = 1/(tfa*(alpm(v)+betm(v)))	
	minf = alpm(v)/(alpm(v)+betm(v))	
}