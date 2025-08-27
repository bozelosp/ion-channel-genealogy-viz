NEURON {
	SUFFIX cat
	USEION ca READ cai, cao
	USEION Ca WRITE iCa VALENCE 2
	RANGE gcatbar, iCa
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	tBase = 23.5  (degC)
	gcatbar = 0 (S/cm2) 
	ki = 0.001 (mM)
	tfa = 1 
	tfi = 0.68 
}

ASSIGNED {
	v (mV)
	celsius (degC)	
	cai (mM)		
	cao (mM)		
	eca (mV)		
	iCa (mA/cm2)
	gcat (S/cm2) 
	minf (1)
	hinf (1)
	taum (ms)
	tauh (ms)
}

STATE {	m h }

BREAKPOINT {
	SOLVE states METHOD cnexp
	gcat = gcatbar*pow(m, 2)*h*h2(cai)	
	iCa = gcat*ghk(v, cai, cao)		
}

DERIVATIVE states { 
	rates(v)
	m' = (minf - m)/taum
	h' = (hinf - h)/tauh
}

INITIAL {
	rates(v)
	m = minf
	h = hinf
}

FUNCTION h2(cai(mM)) {
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
	}else{
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

FUNCTION alph(v (mV)) (/ms) {
	alph = 1.6e-4(/ms)*exp(-(v+57(mV))/19(mV))
}

FUNCTION beth(v (mV)) (/ms) {
	beth = 1(/ms)/(exp(-(v-15(mV))/10(mV)) + 1.0)
}

FUNCTION alpm(v (mV)) (/ms) {
	alpm = 0.1967(/ms)*vtrap(-(v-19.88(mV)), 10.0(mV))
}

FUNCTION betm(v (mV)) (/ms) {
	betm = 0.046(/ms)*exp(-v/22.73(mV))
}

PROCEDURE rates(v (mV)) { 

	taum = 1/(tfa*(alpm(v) + betm(v)))	
	minf = alpm(v)/(alpm(v)+betm(v))	

	tauh = 1/(tfi*(alph(v) + beth(v)))	
	hinf = alph(v)/(alph(v)+beth(v))	
}