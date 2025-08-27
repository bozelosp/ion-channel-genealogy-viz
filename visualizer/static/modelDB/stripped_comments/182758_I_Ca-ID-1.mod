NEURON {
	SUFFIX Ca
	NONSPECIFIC_CURRENT I
	RANGE g0, v0, theta_s, sigma_s
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	v 		(mV)
	g0 		(S/cm2)
	v0 		(mV)

	theta_s
	sigma_s
}

ASSIGNED {
	I 		(mA/cm2)
	s_inf
}

INITIAL {
	rates(v)
}

BREAKPOINT {
	rates(v)

 	I=g0*s_inf*s_inf*(v-v0)
}

PROCEDURE rates(v(mV)) {  
	

UNITSOFF
	s_inf = 1/(1+exp(-(v-theta_s)/sigma_s))
}

UNITSON