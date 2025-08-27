NEURON {
	POINT_PROCESS syn
	POINTER vg
	NONSPECIFIC_CURRENT I
	RANGE g0, v0, theta_g, theta_Hg, sigma_Hg, alpha, beta, s, H_inf,vh
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
	vg		(mV)

	theta_g
	theta_Hg
	sigma_Hg
	alpha 	(ms-1)
	beta 	(ms-1)
}

ASSIGNED {
	I 		(mA/cm2)
	vh 		(mV)
	H_inf
}

STATE {
	s
}

INITIAL {
	rates(v)
	s = 0
}

BREAKPOINT {
 	SOLVE states METHOD cnexp

 	I=g0*(v-v0)*s
}

DERIVATIVE states {
	rates(v)
	s' = alpha*H_inf*(1-s)-beta*s
}

PROCEDURE rates(v(mV)) {  
	

UNITSOFF
	vh = vg-theta_g
	H_inf = 1/(1+exp(-(vh-theta_Hg)/sigma_Hg))
}

UNITSON