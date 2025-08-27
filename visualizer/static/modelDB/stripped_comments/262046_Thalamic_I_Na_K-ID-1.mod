UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX thalamic_i_na_k
	USEION na WRITE ina				
	USEION k WRITE ik
	RANGE g_Na, i_na				
	RANGE g_K, i_k					
}

PARAMETER {
	ena = 50 (mV)
	ek = -90 (mV)
	i_na = 0.0 (mA/cm2)
	g_Na = 0.3 (S/cm2)
	i_k = 0.0 (mA/cm2)
	g_K = 0.5 (S/cm2)
}

ASSIGNED {
	v (mV)
	ina (mA/cm2)
	ik (mA/cm2)
	h_inf
	tau_h (ms)
	m_inf
}

STATE {
	h 
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = g_Na*m_inf*m_inf*m_inf*h*(v - ena)
	ik = g_K*(0.75*(1-h))*(0.75*(1-h))*(0.75*(1-h))*(0.75*(1-h))*(v - ek)
	i_na = ina 							
	i_k = ik							
}

UNITSOFF

INITIAL {
	settables(v)
	h = h_inf
}

DERIVATIVE states {
	settables(v)
	h' = (h_inf - h)/tau_h
}

PROCEDURE settables(v) {
	TABLE h_inf, m_inf, tau_h FROM -100 TO 100 WITH 400
	
	h_inf = 1/(1+exp((v+41)/4))
	tau_h = 1/((0.128*exp(-(v+46)/18))+(4/(1+exp(-(v+23)/5))))
	m_inf = 1/(1+exp(-(v+37)/7))
}

UNITSON