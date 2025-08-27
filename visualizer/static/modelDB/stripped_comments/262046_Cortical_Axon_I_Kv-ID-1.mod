UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX cortical_axon_i_kv
	USEION k WRITE ik				
	RANGE g_Kv, i_Kv				
}

PARAMETER {
	ek = -90 (mV)
	i_Kv = 0.0 (mA/cm2)				
	g_Kv = 0.6e-3 (S/cm2)
	Q_s = 3.209						
}

ASSIGNED {
	v (mV)
	ik (mA/cm2)
	alpha_n
	beta_n
	n_inf
	tau_n (ms)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = g_Kv*n*(v - ek)
	i_Kv = ik 						
}

UNITSOFF

INITIAL {
	settables(v)
	n = n_inf
}

DERIVATIVE states {
	settables(v)
	n' = (n_inf-n)/tau_n
}

PROCEDURE settables(v) {
	TABLE alpha_n, beta_n, n_inf, tau_n FROM -100 TO 100 WITH 400
	
	alpha_n = Q_s*0.01*vtrap(-(v-30),9)
	beta_n = Q_s*0.002*vtrap((v-30),9)
	
	n_inf = alpha_n/(alpha_n+beta_n)
	tau_n = 1/(alpha_n+beta_n)
}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y)-1)
	}
}

UNITSON