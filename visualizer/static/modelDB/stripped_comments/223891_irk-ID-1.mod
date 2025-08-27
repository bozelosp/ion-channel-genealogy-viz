INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX irk
	USEION k READ ek WRITE ik 
  RANGE gkbar, m_inf, tau_m 
  GLOBAL shift
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v		(mV)
	ek		(mV)
	gkbar= 0.0005	(mho/cm2)
  tau_m = 1.0
  shift = 0
}

STATE {
	m
}

ASSIGNED {
	ik		(mA/cm2)
	m_inf
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
 	ik = gkbar * m^3 * (v-ek)
}

DERIVATIVE states { 
	evaluate_fct(v)
	m'= (m_inf-m) / tau_m
}

INITIAL {
	evaluate_fct(v)
	m = m_inf
}

PROCEDURE evaluate_fct(v(mV)) {

	m_inf = 1/(1 + exp( (v - ek - 15 + shift)/10) )
}