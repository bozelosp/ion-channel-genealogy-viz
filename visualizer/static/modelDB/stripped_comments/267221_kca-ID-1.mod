NEURON {
	SUFFIX kca
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gbar, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	gbar = 0.01 (S/cm2)
	beta = 0.03 (1/ms)	
	cac = 0.025 (mM)	
	taumin = 0.5 (ms)	
	q10 = 3 (1) 		
}

ASSIGNED {
	v (mV)
	celsius (degC)
	ek (mV)
	cai (mM)	
	gk (S/cm2)
	ik (mA/cm2)
	minf
	taum (ms)
}

STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	gk = gbar*pow(m, 3)		
	ik = gk*(v - ek)	
}

DERIVATIVE states { 
	rates(cai)
	m' = (minf - m) / taum
}

INITIAL {
	rates(cai)
	m = minf
}

PROCEDURE rates(cai (mM)) {
	LOCAL car, tadj
	tadj = q10^((celsius-22.0(degC))/10(degC))		
	car = (cai/cac)^2
	minf = car / ( 1 + car )			
	taum =  1 / beta / (1 + car) / tadj	
	if (taum < taumin) {taum = taumin}
}