UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX cortical_soma_i_m
	USEION k WRITE ik				
	RANGE g_M, i_M					
}

PARAMETER {
	ek = -100 (mV)
	i_M = 0.0 (mA/cm2)				
	g_M = 7e-5 (S/cm2)
	tau_max = 1000 (ms)
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	p_inf
	tau_p (ms)
}

STATE {
	p
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = g_M*p*(v - ek)
	i_M = ik 						
}

UNITSOFF

INITIAL {
	settables(v)					
	p = p_inf
}

DERIVATIVE states {
	settables(v)
	p' = (p_inf - p)/tau_p
}

PROCEDURE settables(v) {
	TABLE p_inf, tau_p DEPEND tau_max FROM -100 TO 100 WITH 400
	
	p_inf = 1/(1+exp(-(v+35)/10))
	tau_p = tau_max/(3.3*exp((v+35)/20)+exp(-(v+35)/20))
}

UNITSON