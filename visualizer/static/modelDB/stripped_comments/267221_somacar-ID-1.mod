NEURON {
	SUFFIX somacar
	USEION ca READ eca WRITE ica
	RANGE gcabar, ica
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER { 
	gcabar = 0 (S/cm2) 
}

STATE { 
	m
	h
}

ASSIGNED { 
	v (mV)
	eca (mV)	
	ica (mA/cm2)
	minf (1)
	hinf (1)
	taum (ms)
	tauh (ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gcabar*pow(m, 3)*h*(v - eca)
}

DERIVATIVE states {
	rates(v)
	m' = (minf - m)/taum
	h' = (hinf - h)/tauh
}

INITIAL {
	rates(v)
	m = 0	
	h = 1	
}

PROCEDURE rates(v (mV)) {
	minf = 1 / (1 + exp((v + 55(mV))/(-3(mV))))	
	hinf = 1 / (1 + exp((v + 57(mV))/(1(mV))))	
	taum = 100.0
	tauh = 5.0
}