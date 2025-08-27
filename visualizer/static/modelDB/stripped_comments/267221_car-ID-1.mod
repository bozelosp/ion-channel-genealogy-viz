NEURON {
	SUFFIX car
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

STATE {	m h } 

ASSIGNED { 
	v (mV)
	celsius (degC)
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
	m = minf	
	h = hinf	
}

PROCEDURE rates(v (mV)) {
	minf = 1 / (1 + exp(-(v+48.5(mV))/3(mV)))	
	taum = 50

	hinf = 1 / (1 + exp((v+53(mV))/(1(mV))))	
	tauh = 5
}