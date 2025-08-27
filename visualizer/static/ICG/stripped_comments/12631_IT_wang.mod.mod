INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX it
	USEION ca READ cai,cao WRITE ica
        RANGE gcabar, m_inf, tau_m, alph1, alph2, KK, shift, carev

}

UNITS {
	(molar)	= (1/liter)
	(mM)	= (millimolar)
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)
	celsius = 36	(degC)

	gcabar	= .0008	(mho/cm2)
	shift	= 2	(mV)
	cai	= 2.4e-4 (mM)		
	cao	= 2	(mM)
}

STATE {
	m h d
}

ASSIGNED {
	ica 	(mA/cm2)
        carev 	(mV)
	m_inf
	tau_m	(ms)
	alph1	(/ms)
	alph2	(/ms)
	KK
	phi_m
	phi_h
}

BREAKPOINT {
	SOLVE states METHOD adams
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar * m*m*m*h * (v - carev)
}

DERIVATIVE states { 
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = alph1 * ( (1-h-d) - KK * h )
        d' = alph2 * ( KK * (1-h-d) - d )
}

UNITSOFF
INITIAL {
	evaluate_fct(v)
	m = m_inf
        h = 1./(1. + KK + KK*KK)
	d = KK*KK/(1 + KK + KK*KK)




	phi_m = 5.0 ^ ((celsius-24)/10)
	phi_h = 3.0 ^ ((celsius-24)/10)
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL tau2








	m_inf = 1 / ( 1 + exp(-(v+shift+63)/7.8) )
	tau_m = m_inf * ( 1.7 + exp(-(v+shift+28.8)/13.5) ) / phi_m
        alph1 = phi_h * exp(-(v+shift+160.3)/17.8)
        KK    = sqrt( 0.25 + exp((v+shift+83.5)/6.3) ) - 0.5
	tau2  = 240.0 / ( 1 + exp((v+shift+37.4)/30) ) / phi_h
        alph2 = 1 / ( tau2 * (KK+1) )
}
UNITSON