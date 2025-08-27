UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

NEURON {
	SUFFIX IA
	USEION k READ ek WRITE ik
	RANGE gkAbar, ik
}

PARAMETER {
	gkAbar = 0.0165 (S/cm2)	
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	ek (mV)
	ainf
	binf
	tau_b (ms)
	tau_a (ms)
}

STATE {
	a
	b
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkAbar*a*b*(v - ek)
}

DERIVATIVE states { 
	
	
	rates(v)
	a' = (ainf - a)/(tau_a)
	b' = (binf - b)/(tau_b)
}

INITIAL {
	rates(v)
	a = ainf
	b = binf
}

PROCEDURE rates(v (mV)) {
	
	
	LOCAL alpha_b, beta_b

	alpha_b = 0.000009(/ms) / exp((v - 26(mV))/18.5(mV))
	beta_b = 0.014 (/ms) / (exp(-(v + 70(mV))/11(mV)) + 0.2)
	tau_b = 1/(alpha_b + beta_b)
	binf = 1/(1 + exp((v + 71(mV))/7.3(mV)))

	tau_a = 5
	ainf = 1/(1 + exp(-(v + 14(mV))/20.6(mV)))
}