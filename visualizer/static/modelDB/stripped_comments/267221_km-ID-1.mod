NEURON {
	SUFFIX km
	USEION k READ ek WRITE ik
	RANGE gbar, ik
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
}

PARAMETER {
	gbar = 0.03 (S/cm2)
	tha  = -30 (mV)			
	qa   = 9 (mV)			
	Ra   = 0.001 (/ms/mV)	
	Rb   = 0.001 (/ms/mV)	
	temp = 23	(degC)		
	q10  = 2.3				
}

ASSIGNED {
	v (mV)
	gk (S/cm2)
	ek (mV)
	celsius (degC)
	ik (mA/cm2)
	tadj
	ninf
	taun (ms)
}

STATE {
	n
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	tadj = q10^((celsius - temp)/10(degC))  
	gk = tadj*gbar*n
	ik = gk * (v - ek)
}

DERIVATIVE states { 
	rates(v)
	n' = (ninf - n)/taun
}

INITIAL { 
	rates(v)
	n = 0
}

FUNCTION vtrap(x (mV), y (mV)) (1) {
	
	if (fabs(x/y) < 1e-6) {
		vtrap = 1(/mV)*y*(1 - x/y/2)
	} else {  
		vtrap = 1(/mV)*x/(exp(x/y) - 1)
	}
}

FUNCTION alpha(v (mV)) (/ms) { 
	alpha = Ra * vtrap(-(v - tha), qa)
}

FUNCTION beta(v (mV)) (/ms) { 
	beta = Rb * vtrap(v - tha, qa)
}

PROCEDURE rates(v (mV)) {
	
	
	taun = 1/(alpha(v) + beta(v))
	ninf = alpha(v)*taun
}