UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX interneuron_i_k
	USEION k WRITE ik				
	RANGE g_K, i_K					
}

PARAMETER {
	ek = -100 (mV)
	i_K = 0.0 (mA/cm2)				
	g_K = 0.01 (S/cm2)
	V_T = -55(mV)
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	alpha_n
	beta_n
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = g_K*n*n*n*n*(v - ek)
	i_K = ik 						
}

UNITSOFF

INITIAL {
	settables(v)					
	n = 0
}

DERIVATIVE states {
	settables(v)
	n' = alpha_n*(1-n)-beta_n*n
}

PROCEDURE settables(v) {
	TABLE alpha_n, beta_n DEPEND V_T FROM -100 TO 100 WITH 400
	alpha_n = 0.032 * vtrap(-(v-V_T-15), 5)
	beta_n = 0.5*exp(-(v-V_T-10)/40)
}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y)-1)
	}
}


UNITSON