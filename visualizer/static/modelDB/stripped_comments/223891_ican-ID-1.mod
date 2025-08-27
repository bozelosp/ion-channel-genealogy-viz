INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ican
	USEION n READ en WRITE in VALENCE 1
	USEION ca READ cai
        RANGE gbar, m_inf, tau_m
	GLOBAL beta, cac, taumin
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


PARAMETER {
	v		(mV)
	celsius	= 36	(degC)
	en	= -20	(mV)		
	cai 	= 2.4e-4 (mM)		
	gbar	= 0.00025 (mho/cm2)
	beta	= 0.002	(1/ms)		
	cac	= 0.01	(mM)		
	taumin	= 0.1	(ms)		
}


STATE {
	m
}

ASSIGNED {
	in	(mA/cm2)
	m_inf
	tau_m	(ms)
	tadj
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	in = gbar * m*m * (v - en)
}

DERIVATIVE states { 
	evaluate_fct(v,cai)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {




	tadj = 3.0 ^ ((celsius-22.0)/10)

	evaluate_fct(v,cai)
	m = m_inf
}


PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL alpha2

	alpha2 = beta * (cai/cac)^2

	tau_m = 1 / (alpha2 + beta) / tadj
	m_inf = alpha2 / (alpha2 + beta)

        if(tau_m < taumin) { tau_m = taumin } 	
}
UNITSON