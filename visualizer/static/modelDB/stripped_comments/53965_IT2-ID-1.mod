INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX it2
	USEION ca READ cai, cao WRITE ica
	RANGE gcabar, m_inf, tau_m, h_inf, tau_h, shift
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

	gcabar	= .0008	(mho/cm2)
	shift	= 0 	(mV)
	cai	= 2.4e-4 (mM)		
	cao	= 2	(mM)
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






	phi_m = 5.0 ^ (12/10)
	phi_h = 3.0 ^ (12/10)

	evaluate_fct(v)

	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 




	m_inf = 1.0 / ( 1 + exp(-(v+shift+52.5)/7.85) )
	h_inf = 1.0 / ( 1 + exp((v+shift+83.24)/8.51) )

	tau_m = ( 1.97 + 1.0 / ( exp((v+shift+26.5)/15.3) + exp(-(v+shift+147)/27) ) ) / phi_m
	tau_h = 16 + (1942 + exp((v+shift+161)/8.7)) / (1 + exp((v+shift+89.6)/3.7) ) / phi_h
}
UNITSON