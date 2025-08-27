NEURON {
	SUFFIX Ks
	USEION k READ ek WRITE ik
	RANGE gkbar, vtraub, rate_change
	RANGE n_inf
	RANGE tau_n
	RANGE ik 
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	gkbar		= .0003 	(mho/cm2)
	ek				(mV)
	celsius			(degC)
	v               		(mV)
	vtraub	= -55		(mV)	
	rate_change = 0.1			
}

STATE {
	n
}

ASSIGNED {
	ik	(mA/cm2)
	n_inf
	tau_n (ms)
	tadj
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ik  = gkbar * n * (v - ek)
}

DERIVATIVE states {
	evaluate_fct(v)
	n' = (n_inf-n)/tau_n
}

UNITSOFF
INITIAL {



	tadj = 3.0 ^ ((celsius-36)/ 10 )
	evaluate_fct(v)
	n= n_inf
}

PROCEDURE evaluate_fct(v(mV)) { LOCAL a,b,v2

	v2 = v - vtraub 

	a = 0.032 * (15-v2) / ( exp((15-v2)/5) - 1)
	b = 0.5 * exp((10-v2)/40)
	tau_n = 1 / (a + b) / (tadj*rate_change)
	n_inf = a / (a + b)
}

UNITSON