NEURON {
	SUFFIX calH
	USEION ca READ eca WRITE ica
	RANGE gcalbar, ica
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {	
	gcalbar = 0 (S/cm2)	
}

ASSIGNED {	
	v (mV)
	celsius (degC)
	ica (mA/cm2)
	eca (mV)	
	minf (1)
	hinf (1)
	taum (ms)
	tauh (ms)
}

STATE {	
	m
	h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = gcalbar*pow(m, 3)*h*(v - eca)
}

DERIVATIVE states {
	rates (v)
	m' = (minf-m)/taum
	h' = (hinf-h)/tauh
}

INITIAL {
	rates(v)
	m = minf	
	h = hinf	
}

PROCEDURE rates(v (mV)) {
	minf = 1 / (1 + exp(-(v+37(mV))/1(mV)))		
	taum = 3.6

	hinf = 1 / (1 + exp((v+41(mV))/0.5(mV)))	
	tauh = 29
}