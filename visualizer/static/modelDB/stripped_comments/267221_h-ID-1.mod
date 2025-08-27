NEURON {
	SUFFIX h
	NONSPECIFIC_CURRENT ih
	RANGE  ghbar, vhalf, eh, ih
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	ghbar = 0.0 (S/cm2)  
	eh = -10 (mV)
	K = 3.5 (mV)
	vhalf = -90 (mV)    
}

ASSIGNED {
	v (mV)
	ih (mA/cm2)
	ninf
	taun (ms)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ih = ghbar*n*(v-eh)
}

DERIVATIVE states {
	rates(v)
	n' = (ninf - n)/taun
}

INITIAL {
	rates(v)
	n = ninf
}

PROCEDURE rates(v (mV)) {
	if (v > -30) {
		taun = 1
	} else {
		taun = 2(ms)*(1/(exp(-(v+145(mV))/17.5(mV)) + exp((v+16.8(mV))/16.5(mV))) + 10) 
	}
	ninf = (1 / (1 + exp((v - vhalf)/K)))	
}