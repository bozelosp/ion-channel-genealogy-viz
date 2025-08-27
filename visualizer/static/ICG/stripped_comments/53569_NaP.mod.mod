NEURON {
	SUFFIX NaP
	USEION na READ ena WRITE ina
	RANGE gnabar, vtraub, vsm, vsh, gamma
	RANGE m_inf, h_inf
	RANGE tau_m, tau_h
	RANGE ina 
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gnabar	= .00029 	(mho/cm2)
	ena				(mV)
	celsius			(degC)
	v               		(mV)
	vtraub	= -55		(mV)	
	vsm		= -2 (mV)	
	vsh		= -5 (mV)	
	gamma		= 0.5		
}

STATE {
	m h
}

ASSIGNED {
	ina	(mA/cm2)
	m_inf
	h_inf
	tau_m (ms)
	tau_h (ms)
	tadj
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar * m*h * (v - ena)
}

DERIVATIVE states {
	evaluate_fct(v)
	m' = (m_inf-m)/tau_m
	h' = (h_inf-h)/tau_h
}

UNITSOFF
INITIAL {



	tadj = 3.0 ^ ((celsius-36)/ 10 )
	evaluate_fct(v)
	m= m_inf
	h= h_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v - vtraub 

	a = 0.32 * (vsm+13-v2) / ( exp((vsm+13-v2)/4) - 1)
	b = 0.28 * (vsm+v2-40) / ( exp((vsm+v2-40)/5) - 1)
	tau_m = 1 / (a + b) / tadj
	m_inf = a / (a + b)

	a = 0.128 * exp((vsh+17-v2)/18)
	b = 4 / ( 1 + exp((vsh+40-v2)/5*gamma))
	tau_h = 1 / (a + b) / tadj
	h_inf = a / (a + b)
}

UNITSON