INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX im
	USEION k READ ek WRITE ik
        RANGE gkbar, m_inf, tau_m, i
	GLOBAL taumax

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}


PARAMETER {
	v		(mV)
	celsius = 36    (degC)
	ek		(mV)
	gkbar	= 1e-6	(mho/cm2)
	taumax	= 1000	(ms)		
}



STATE {
	m
}

ASSIGNED {
	ik	(mA/cm2)
	i	(mA/cm2)
	m_inf
	tau_m	(ms)
	tau_peak	(ms)
	tadj
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar * m * (v - ek)
	i = ik
}

DERIVATIVE states { 
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {
	evaluate_fct(v)
	m = 0



        tadj = 2.3 ^ ((celsius-36)/10)
	tau_peak = taumax / tadj
}

PROCEDURE evaluate_fct(v(mV)) {

	m_inf = 1 / ( 1 + exptable(-(v+35)/10) )
	tau_m = tau_peak / ( 3.3 * exptable((v+35)/20) + exptable(-(v+35)/20) )
}
UNITSON


FUNCTION exptable(x) { 
	TABLE  FROM -25 TO 25 WITH 10000

	if ((x > -25) && (x < 25)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}