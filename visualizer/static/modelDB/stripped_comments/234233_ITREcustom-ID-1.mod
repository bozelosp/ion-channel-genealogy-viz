INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX itrecustom
	USEION ca READ cai, cao WRITE ica
	RANGE gcabar, m_inf, tau_m, h_inf, tau_h, shift, qm, qh, taubase
}

UNITS {
	(molar) = (1/liter)
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)

	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}

PARAMETER {
	v		(mV)
	celsius	= 36	(degC)

	gcabar	= .003	(mho/cm2)
	shift	= 0 	(mV)
	taubase (mV)
	cai	= 2.4e-4 (mM)		
	cao	= 2	(mM)
	qm	= 2.5
	qh 	= 2.5
}

STATE {
	m h
}

ASSIGNED {
	ica	(mA/cm2)
	carev	(mV)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	phi_m
	phi_h
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar * m*m*h * (v-carev)     

}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {





	phi_m = qm ^ ((celsius-24)/10)
	phi_h = qh ^ ((celsius-24)/10)

	evaluate_fct(v)
	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 






	m_inf = 1.0 / ( 1 + exp(-(v+shift+50)/7.4) )
	h_inf = 1.0 / ( 1 + exp((v+shift+78)/5.0) )

	tau_m = ( 3 + 1.0 / ( exp((v+shift+25)/10) + exp(-(v+shift+100)/15) ) )
	tau_h = ( taubase + 1.0 / ( exp((v+shift+46)/4) + exp(-(v+shift+405)/50) ))


	tau_h = tau_h / phi_h
	tau_m = tau_m / phi_m
}
UNITSON