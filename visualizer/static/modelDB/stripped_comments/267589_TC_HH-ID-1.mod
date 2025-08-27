INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX TC_HH
	USEION na READ ena WRITE ina
	USEION k READ ek WRITE ik
	RANGE gna_max, gk_max, vtraub, vtraub2, i_rec
	RANGE m_inf, h_inf, n_inf
	RANGE tau_m, tau_h, tau_n
	RANGE m_exp, h_exp, n_exp
	RANGE ina, ik
	RANGE m, h, n
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S)  = (siemens)
}

PARAMETER {
	gna_max	= 1.0e-1 	(S/cm2) 
	gk_max	= 1.0e-1 	(S/cm2) 

	celsius         (degC)
	dt              (ms)
	v               (mV)
	vtraub = -55.5   
	vtraub2 = -45.5  
}

STATE {
	m h n
}

ASSIGNED {
	ina	(mA/cm2)
	ik	(mA/cm2)
	ena	(mV)
	ek	(mV)
	i_rec	(mA/cm2)
	m_inf
	h_inf
	n_inf
	tau_m
	tau_h
	tau_n
	m_exp
	h_exp
	n_exp
	tcorr
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ina   = gna_max * m*m*m*h * (v - ena)
	ik    = gk_max * n*n*n*n * (v - ek)
	i_rec = ina + ik
}


DERIVATIVE states {   
	evaluate_fct(v)
	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
	n' = (n_inf - n) / tau_n
}









UNITSOFF
INITIAL {
	m = 0
	h = 0
	n = 0





	tcorr = 3.0 ^ ((celsius-36)/ 10 )



}


PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2, v3

	v2 = v - vtraub 
	v3 = v - vtraub2 

	if(v2 == 13 || v2 == 40 || v2 == 15 ){
	v2 = v2+0.0001
	}


	
	
	a = 0.32 * vtrap(v2-13,4)
	b = 0.28 * vtrap(40-v2,5)
	tau_m = 1 / (a + b) / tcorr
	m_inf = a / (a + b)

	a = 0.128 * exp((17-v2)/18)
	b = 4 / ( 1 + exp((40-v2)/5) )
	tau_h = 1 / (a + b) / tcorr
	h_inf = a / (a + b)

	
	a = 0.032 * vtrap(v3-15,5)
	b = 0.5 * exp((10-v3)/40)
	tau_n = 1 / (a + b) / tcorr
	n_inf = a / (a + b)

	m_exp = 1 - exp(-dt/tau_m)
	h_exp = 1 - exp(-dt/tau_h)
	n_exp = 1 - exp(-dt/tau_n)

}

FUNCTION vtrap(x,c) { 
	
	if (fabs(x/c) < 1e-6) {
		vtrap = c + x/2
	} else {
		vtrap = x / (1-exp(-x/c)) 
	}
}
UNITSON