INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ican
	USEION other WRITE iother VALENCE 1
	USEION ca READ cai
        RANGE gbar, i
	GLOBAL m_inf, tau_m, beta, cac, taumin, erev
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
	erev = -20	(mV)
	cai 	= .00005	(mM)	
	gbar	= 1e-5	(mho/cm2)
	beta	= 2.5	(1/ms)		
	cac	= 1e-4	(mM)		
	taumin	= 0.1	(ms)		
}


STATE {
	m
}

INITIAL {
	evaluate_fct(v,cai)
	m = m_inf




	tadj = 3.0 ^ ((celsius-22.0)/10)

}

ASSIGNED {
	i	(mA/cm2)
	iother	(mA/cm2)
	m_inf
	tau_m	(ms)
	tadj
}

BREAKPOINT { 
	SOLVE states METHOD runge
	i = gbar * m*m * (v - erev)
	iother = i
}

DERIVATIVE states { 
	evaluate_fct(v,cai)

	m' = (m_inf - m) / tau_m
}

UNITSOFF

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL alpha2

	alpha2 = beta * (cai/cac)^2

	tau_m = 1 / (alpha2 + beta) / tadj
	m_inf = alpha2 / (alpha2 + beta)

        if(tau_m < taumin) { tau_m = taumin } 	
}
UNITSON