INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iC
	USEION k READ ek WRITE ik
	USEION ca READ cai
        RANGE gkbar, m_inf, tau_m, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
}

PARAMETER {
	v		(mV)
	celsius		(degC)
        dt              (ms)
	ek		(mV)
	cai		(mM)
	gkbar= 0.00345	(mho/cm2) 
}

STATE {
	m
}

ASSIGNED {
	ik		(mA/cm2)
	tau_m		(ms)
	m_inf
	tadj
}

BREAKPOINT { 
	SOLVE states 
	ik = gkbar * m * (v - ek)
}







PROCEDURE states() { 
	rates(v,cai)

	m= m + (1-exp(-dt/tau_m))*(m_inf-m)
}

UNITSOFF
INITIAL {
	tadj = 3^((celsius-23.5)/10)
	rates(v,cai)
	m = m_inf
}

PROCEDURE rates( v(mV), cai(mM)) {  LOCAL a,b
	a = 250 * cai * exp(v/24)
	b = 0.1 * exp(-v/24)
	tau_m = 1/(a+b) / tadj
	m_inf = a/(a+b)
}
UNITSON