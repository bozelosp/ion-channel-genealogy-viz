INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iAHP
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
	dt		(ms)
	ek		(mV)
	cai		(mM)
	gkbar= 0.000807	(mho/cm2)

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
	ik = gkbar * m*m * (v-ek)
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
	a = 1200 * cai^2
	b = 0.001
	tau_m = 1/(a+b) / tadj
        m_inf = a/(a+b)
}
UNITSON