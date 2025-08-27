INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iao
	
	USEION k READ ek WRITE ik
        RANGE i, gabar
	GLOBAL m_inf, tau_m, h_inf, tau_h, n_inf, tau_n, shm, shh, shn, hx, nx
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
}

PARAMETER {
	v		(mV)
	celsius	= 36	(degC)
	
	gabar	= 1.0	(mho/cm2)
	shm	= 15 	(mV)
	shh	= 15 	(mV)
	shn	= 15 	(mV)
	hx	= 0 	(mV)
	nx	= 0 	(mV)
}

STATE {
	m h n
}

ASSIGNED {
        ek      (mV)
	ik	(mA/cm2)
	i	(mA/cm2)
	m_inf
	tau_m	(ms)
	h_inf
	tau_h	(ms)
	n_inf
        tau_n	(ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	i = gabar * (m*m*m*m*h * (v-ek) * 0.6 + m*m*m*m*n * (v-ek) * 0.4)
        ik = i
}

DERIVATIVE states {
	evaluate_fct(v)

	m' = (m_inf - m) / tau_m
	h' = (h_inf - h) / tau_h 
	n' = (h_inf - n) / tau_n 
}

UNITSOFF
INITIAL {
	evaluate_fct(v)
	m = m_inf
	h = h_inf
	n = n_inf




}

PROCEDURE evaluate_fct(v(mV)) { 





	m_inf = 1.0 / ( 1 + exp(-(v+shm+34)/7.4) )
	h_inf = 1.0 / ( 1 + exp((v+shh+78)/5.0) )
	n_inf = 1.0 / ( 1 + exp((v+shn+78)/5.0) )

	tau_m = ( 1.0 / ( exp((v+shm+35.8)/19.7) + exp(-(v+shm+79.7)/12.7) ) + 0.37 )  

        if (v < -80-shh) {
		tau_h = ( 1.0 / ( exp((v+shh+46)/5) + exp(-(v+shh+238)/37.5) ) )
	} else {	
	        tau_h = 70 + hx
	}
        if (v < -73-shn) {
                tau_n = ( 1.0 / ( exp((v+shn+46)/5) + exp(-(v+shn+238)/37.5) ) ) 
        } else {
                tau_n = 60 + nx
        }
}
UNITSON