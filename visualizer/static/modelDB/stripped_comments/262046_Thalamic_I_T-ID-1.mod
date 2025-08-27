UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX thalamic_i_t
	USEION ca WRITE ica				
	RANGE g_T, i_t					
}

PARAMETER {
	eca = 0 (mV)
	i_t = 0.0 (mA/cm2)
	g_T = 0.5 (S/cm2)
}

ASSIGNED {
	v (mV)
	ica (mA/cm2)
	r_inf
	tau_r (ms)
	p_inf
}

STATE {
	r 
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ica = g_T*p_inf*p_inf*r*(v - eca)
	i_t = ica 							
}

UNITSOFF

INITIAL {
	settables(v)
	r = r_inf
}

DERIVATIVE states {
	settables(v)
	r' = (r_inf - r)/tau_r
}

PROCEDURE settables(v) {
	TABLE r_inf, p_inf, tau_r FROM -100 TO 100 WITH 400
	
	r_inf = 1/(1+exp((v+84)/4))
	tau_r = (28+exp(-(v+25)/10.5))
	p_inf = 1/(1+exp(-(v+60)/6.2))
}

UNITSON