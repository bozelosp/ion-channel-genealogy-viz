NEURON {
	SUFFIX kca
	USEION k READ ek WRITE ik
	
	USEION ca READ cai
	RANGE  gbar, po, ik
	GLOBAL m_inf, tau_m
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

ASSIGNED {       
	v               (mV)
	celsius         (degC)
	ek              (mV)
	cai             (mM)           
	ik              (mA/cm2)
	po
	m_inf
	tau_m           (ms)
}

PARAMETER {
	gbar    = 0.01   (mho/cm2)
 	taumin  = 100    (ms)            

	b = 0.005 	(mM)

}


STATE {
	m
}

BREAKPOINT { 
	SOLVE states METHOD cnexp

	po = m*m

	ik = gbar*po*(v - ek)    
}

DERIVATIVE states {
	rates(cai)
	m' = (m_inf - m) / tau_m
}


INITIAL {
	rates(cai)
	m = m_inf
}





PROCEDURE rates(cai(mM)) {  LOCAL a 

	a = cai/b
 	m_inf = a/(a+1)
	tau_m = taumin+ 1(ms)*1(mM)/(cai+b)
}