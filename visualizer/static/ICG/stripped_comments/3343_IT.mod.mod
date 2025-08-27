INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX it
	USEION ca READ cai,cao WRITE ica
	GLOBAL q10
	RANGE gcabar, m_inf, tau_m, h_inf, tau_h, shift, carev
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
	gcabar	= 0.002	(mho/cm2)
	q10	= 3			
	shift	= 2 	(mV)		
	cai	= 2.4e-4 (mM)		
	cao	= 2	(mM)
}

STATE {
	h
}

ASSIGNED {
	ica	(mA/cm2)
	carev	(mV)
	m_inf
	tau_m	(ms)			
	h_inf
	tau_h	(ms)
	phi_h
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar * m_inf * m_inf * h * (v-carev)
}

DERIVATIVE castate {
	evaluate_fct(v)

	h' = (h_inf - h) / tau_h
}


UNITSOFF
INITIAL {




	phi_h = q10 ^ ((celsius-24 (degC) )/10 (degC) )

	h = 0
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL Vm

	Vm = v + shift

	m_inf = 1.0 / ( 1 + exp(-(Vm+57)/6.2) )
	h_inf = 1.0 / ( 1 + exp((Vm+81)/4.0) )







	tau_h = 30.8 + (211.4 + exp((Vm+113.2)/5)) / (1 + exp((Vm+84)/3.2))

	tau_h = tau_h / phi_h

}

UNITSON