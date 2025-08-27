NEURON {
	SUFFIX CaT
	USEION ca READ cai, cao WRITE ica
	RANGE gcabar, m_inf, tau_m, h_inf, tau_h, shift
	RANGE vsm, vsh
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
	gcabar = .00003 	(mho/cm2)
	shift	= 2 		(mV)		
	cai			(mM)		
	cao			(mM)
	vsm = 20		(mV)		
	vsh = 20 		(mV)		
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

 	m'= (m_inf-m) / tau_m
     	h'= (h_inf-h) / tau_h
}

UNITSOFF
INITIAL {






	phi_m = 5.0 ^ ((celsius-24)/10)
	phi_h = 3.0 ^ ((celsius-24)/10)	

	evaluate_fct(v)
	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 




	m_inf = 1.0 / ( 1 + exp(-(v+shift+vsm+50)/(7.4)) )
	h_inf = 1.0 / ( 1 + exp((v+shift+vsh+78)/5.0) )

	tau_m = ( 3+1.0 / ( exp((v+shift+25)/10) + exp(-(v+shift+100)/15) ) ) / phi_m
	tau_h = ( 85+1.0 / ( exp((v+shift+46)/4) + exp(-(v+shift+405)/50) ) ) / phi_h
}
UNITSON