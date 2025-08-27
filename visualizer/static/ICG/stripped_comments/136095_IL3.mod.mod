INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ical3
	
        USEION ca READ cai,cao WRITE ica
        RANGE pcabar, g
	GLOBAL 	m_inf, tau_m, sh1, sh2, rat
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


PARAMETER {
	v		(mV)
	celsius	= 36	(degC)
	eca     = 120		(mV)
	cai 	= .00005	(mM)	
	cao 	= 2		(mM)	
	pcabar	= .003	(mho/cm2)	
	sh1 	= 0		 
	sh2	= 0
	rat	= 1
}


STATE {
	m
}

INITIAL {
	tadj = 3 ^ ((celsius-21.0)/10)
	evaluate_fct(v)
	m = m_inf
}


ASSIGNED {
	ica	(mA/cm2)
	g       (mho/cm2)
	m_inf
	tau_m	(ms)
        tadj
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	g = pcabar * m * m
	ica = g * ghk(v, cai, cao)
}

DERIVATIVE states { 
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
PROCEDURE evaluate_fct(v(mV)) {  LOCAL a,b




	a = 1.6 / (1 + exp(-0.072*(v+sh1+5)) )
	b = 0.02 * (v+sh2-1.31) / ( exp((v+sh2-1.31)/5.36) - 1)

	tau_m = 1.0 / (a + b) / tadj
	m_inf = a / (a + b)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(z)*rat
	eci = ci*efun(-z)
	
	
	ghk = (.001)*2*FARADAY*(eci - eco)
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}
UNITSON