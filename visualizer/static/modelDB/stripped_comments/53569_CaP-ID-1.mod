NEURON {
	SUFFIX CaP
	USEION ca READ cai, cao WRITE ica
	RANGE gcabar, m_inf, tau_m, shift
	RANGE delta
	RANGE ica
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
	v			(mV)
	celsius		(degC)
	gcabar = .00002 	(mho/cm2)
	shift	= 2 		(mV)		
	cai			(mM)		
	cao			(mM)
	delta = 60		(mV)		
}

STATE {
	m 
}

ASSIGNED {
	ica	(mA/cm2)
	carev	(mV)
	m_inf
	tau_m	(ms)
	phi_m
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar * m*m * (v-carev)
}

DERIVATIVE castate {
	evaluate_fct(v)

 	m'= (m_inf-m) / tau_m
}

UNITSOFF
INITIAL {






	phi_m = 5.0 ^ ((celsius-24)/10)
	evaluate_fct(v)
	m = m_inf
}

PROCEDURE evaluate_fct(v(mV)) { 



	m_inf = 1.0/(1+exp(-(v+shift+50)/(7.4)))
	tau_m = ( 3+delta+1.0 /(exp((v+shift+25)/10)+exp(-(v+shift+100)/15)))/phi_m
}
UNITSON