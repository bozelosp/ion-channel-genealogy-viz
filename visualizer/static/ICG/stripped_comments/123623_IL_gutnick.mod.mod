INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ical
	USEION ca READ cai, cao, eca WRITE ica
        RANGE gcabar, alpha_m, beta_m, alpha_h, beta_h, m, h, carev
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
	eca		(mV)
	cai 	= .00024 (mM)		
	cao 	= 2	(mM)		
	gcabar	= 1e-4	(mho/cm2)	
}


STATE {
	m
	h
}

ASSIGNED {
	ica	(mA/cm2)		
	carev	(mV)			
	alpha_m	(/ms)			
	beta_m	(/ms)
	alpha_h	(/ms)
	beta_h	(/ms)
	tadj
}


BREAKPOINT { 
	SOLVE states METHOD cnexp 
	carev = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar * m * m * h * (v-carev)
}

DERIVATIVE states { 
	evaluate_fct(v)

	m' = alpha_m * (1-m) - beta_m * m
	h' = alpha_h * (1-h) - beta_h * h
}


UNITSOFF

INITIAL {
	evaluate_fct(v)



}

PROCEDURE evaluate_fct(v(mV)) {

	

	alpha_m = 0.055 * (-27-v) / (exp((-27-v)/3.8) - 1)
	beta_m = 0.94 * exp((-75-v)/17)

	alpha_h = 0.000457 * exp((-13-v)/50)
	beta_h = 0.0065 / (exp((-15-v)/28) + 1)
}

UNITSON