INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iahp
	USEION k READ ek WRITE ik
	USEION ca READ cai
        RANGE gkbar, i
	GLOBAL beta, cac, m_inf, tau_m
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
	
	
	gkbar	= .001	(mho/cm2)
	beta	= 2.5	(1/ms)		
	cac	= 1e-4	(mM)		
	taumin	= 1	(ms)		
}


STATE {
	m
}

ASSIGNED {
        cai (mM)
        ek (mV)
	ik	(mA/cm2)
	i	(mA/cm2)
	m_inf
	tau_m	(ms)
	tadj
}

BREAKPOINT { 
	SOLVE states METHOD runge
	i = gkbar * m*m * (v - ek)
	ik = i
}

DERIVATIVE states { 
	evaluate_fct(v,cai)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {
	evaluate_fct(v,cai)
	m = m_inf




	tadj = 3 ^ ((celsius-22.0)/10)

}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL car

	car = (cai/cac)^2

	m_inf = car / ( 1 + car )
	tau_m = 1 / beta / (1 + car) / tadj

        if(tau_m < taumin) { tau_m = taumin } 	
}
UNITSON