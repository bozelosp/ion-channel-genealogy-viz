INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iahp
	
	USEION k READ ek WRITE ik
        
	USEION ca READ cai
        RANGE gkbar, i, g, ratc, ratC, minf, taum, Cai
	GLOBAL beta, cac, m_inf, tau_m, x
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
	
	Cai 	= 0.0	(mM)		
	
	gkbar	= .001	(mho/cm2)
	beta	= 2.5	(1/ms)		
	cac	= 1e-4	(mM)		
	taumin	= 1	(ms)		
        ratc    = 1
        ratC    = 1
        x       = 2
}


STATE {
	m
}
ASSIGNED {
        cai (mM)
        ek (mV)
	ik 	(mA/cm2)
	i	(mA/cm2)
	g       (mho/cm2)
	m_inf
	tau_m	(ms)
	minf
        taum
	tadj
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
        minf=m_inf
        taum=tau_m
	g = gkbar * m*m
	i = g * (v - ek)
	ik  = i
}

DERIVATIVE states { 
	evaluate_fct(v,Cai,cai)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {




	

	tadj = 3 ^ ((celsius-22.0)/10)

	evaluate_fct(v,Cai,cai)
	m = m_inf
        minf=m_inf
        taum=tau_m
}

PROCEDURE evaluate_fct(v(mV),Cai(mM), cai(mM)) {  LOCAL car, tcar

        tcar = ratC*Cai + ratc*cai
	car = (tcar/cac)^x

	m_inf = car / ( 1 + car )
	tau_m = 1 / beta / (1 + car) / tadj

        if(tau_m < taumin) { tau_m = taumin } 	
}
UNITSON