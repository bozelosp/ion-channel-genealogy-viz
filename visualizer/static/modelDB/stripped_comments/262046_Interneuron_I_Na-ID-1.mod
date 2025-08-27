UNITS {
 (mV) = (millivolt)
 (mA) = (milliamp)
 (S) = (siemens)
}

NEURON {
	SUFFIX interneuron_i_na
	USEION na WRITE ina				
	RANGE g_Na, i_Na				
}

PARAMETER {
	ena = 50 (mV)
	i_Na = 0.0 (mA/cm2)				
	g_Na = 0.05 (S/cm2)
	V_T = -55(mV)
}

ASSIGNED {
	v (mV)
	ina (mA/cm2)
	alpha_m
	beta_m
	alpha_h
	beta_h
}

STATE {
	m h
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = g_Na*m*m*m*h*(v - ena)
	i_Na = ina 						
}

UNITSOFF

INITIAL {
	settables(v)					
	m = 0
	h = 0
}

DERIVATIVE states {
	settables(v)
	m' = alpha_m*(1-m)-beta_m*m
	h' = alpha_h*(1-h)-beta_h*h
}

PROCEDURE settables(v) {
	TABLE alpha_m, beta_m, alpha_h, beta_h DEPEND V_T FROM -100 TO 100 WITH 400
	
	alpha_m = 0.32*vtrap(-(v-V_T-13),4)
	beta_m = 0.28*vtrap((v-V_T-40),5)
	alpha_h = 0.128*exp(-(v-V_T-17)/18)
	beta_h = 4/(1+exp(-(v-V_T-40)/5))
}

FUNCTION vtrap(x,y) {
	if (fabs(x/y) < 1e-6) {
		vtrap = y*(1 - x/y/2)
	}else{
		vtrap = x/(exp(x/y)-1)
	}
}

UNITSON