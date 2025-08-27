INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX itre
	USEION ca READ eca, cai, cao WRITE ica
	RANGE gmax, m_inf, tau_m, h_inf, tau_h, carev, shift, i
        GLOBAL exptemp, q10m, q10h
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
	gmax	= .003	(mho/cm2)
	shift	= 2 	(mV)
	q10m	= 2.5
	q10h	= 2.5
        exptemp = 24
        cao
        cai

}

STATE {
	m h
}

ASSIGNED {
	i	(mA/cm2)  
	ica	(mA/cm2)
        eca (mV)
	
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	phim
        phih
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	
	i = gmax * m*m*h * (v-eca)
        ica=i
}

DERIVATIVE states {
	mh(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h
}

UNITSOFF
INITIAL {





	phim = q10m ^ ((celsius-exptemp)/10)
	phih = q10h ^ ((celsius-exptemp)/10)

	mh(v)
	m = m_inf
	h = h_inf
}

PROCEDURE mh(v(mV)) { 




	m_inf = 1.0 / ( 1 + exp(-(v+shift+50)/7.4) )
	h_inf = 1.0 / ( 1 + exp((v+shift+78)/5.0) )

	tau_m = ( 1 + 0.33 / ( exp((v+shift+25)/10) + exp(-(v+shift+100)/15) ) ) / phim


	tau_h = ( 85 + 1.0 / ( exp((v+shift+46)/4) + exp(-(v+shift+405)/50) ) ) / phih
}
UNITSON