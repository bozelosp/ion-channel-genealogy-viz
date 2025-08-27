INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iCat1m3h
	USEION ca READ eca, cai, cao WRITE ica
	RANGE gcabar, m_inf, tau_m, h_inf, tau_h, shift, i
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
	celsius	= 35	(degC)
	gcabar	= .0008	(mho/cm2)
	cai	= 8e-5 (mM)		
	cao	= 2	(mM)
}

STATE {
	m h
}

ASSIGNED {
        eca     (mV)
	ica	(mA/cm2)
	
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	i	(mA/cm2)
}

BREAKPOINT {
	SOLVE castate METHOD cnexp
	


	ica = gcabar * m*m*m*h * (v-eca)
	i = ica		
}

DERIVATIVE castate {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {


	evaluate_fct(v)

	m = m_inf
	h = h_inf
}

PROCEDURE evaluate_fct(v(mV)) { 

	m_inf = 1.0/(1+exp((-51.44-v)/7.23))
	h_inf = 1.0 / ( 1 + exp((v+73.43)/6.04) )

	tau_m = 2.29 + 1.0 / ( exp((v+68.03)/-27.68) + exp((v+39.08)/2.74) ) 
	if (v < -50)
	{tau_h = exp((v+770)/162.5)}
	else
	{tau_h = 37 + exp((v+27)/-6)}
}

UNITSON