INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX icanv
	USEION n READ en WRITE in VALENCE 1
	USEION ca READ cai
        RANGE gbar
	GLOBAL 	m_inf, tau_m, cac, taumin, vact, vtau
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
	en		(mV)
	cai 	= .00005	(mM)	
	gbar	= 1e-5	(mho/cm2)
	cac	= 1e-4	(mM)		
	taumin	= 0.1	(ms)		
	vact	= -64	(mV)		
	vtau	= -92	(mV)		
}


STATE {
	m
}

INITIAL {
	evaluate_fct(v,cai)
	m = m_inf
}


ASSIGNED {
	in	(mA/cm2)
	m_inf
	tau_m	(ms)
}

BREAKPOINT { 
	SOLVE states
	in = gbar * m*m * (v - en)
}

DERIVATIVE states { 
	evaluate_fct(v,cai)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL cc,tadj






	tadj = 3 ^ ((celsius-22.0)/10)

	cc = (cai/cac)^2

	m_inf = 1 / (1 + exp(-(v-vact)/2) / cc )

	tau_m = exp((v-vtau)/4) / (1 + cc*exp((v-vact)/2) ) / tadj

        if(tau_m < taumin) { tau_m = taumin } 	
}
UNITSON