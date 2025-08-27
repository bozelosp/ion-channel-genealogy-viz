INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ical
	USEION Ca READ Cai, Cao WRITE iCa VALENCE 2
      RANGE pcabar, g
	GLOBAL 	m_inf, taum, sh1, sh2
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
	eCa     = 120		(mV)
	Cai 	= .00005	(mM)	
	Cao 	= 2		(mM)	
	pcabar	= 9e-4	(mho/cm2)
	sh1 	= -17		 
	sh2	= -7		 
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
	iCa	(mA/cm2)
	g       (mho/cm2)
	m_inf
	taum	(ms)
      tadj
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	g = pcabar * m * m
	iCa = g * ghk(v, Cai, Cao)
}

DERIVATIVE states { 
	evaluate_fct(v)
	m' = (m_inf - m) / taum
}


UNITSOFF
PROCEDURE evaluate_fct(v(mV)) {  LOCAL a,b



	a = 1.6 / (1 + exp(-0.072*(v+sh1+5)) )
	b = 0.02 * (v+sh2-1.31) / ( exp((v+sh2-1.31)/5.36) - 1)
	taum = 1.0 / (a + b) / tadj
	m_inf = a / (a + b)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (.001 coul/cm3) {
	LOCAL z, eci, eco
	z = (1e-3)*2*FARADAY*v/(R*(celsius+273.15))
	eco = co*efun(z)
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