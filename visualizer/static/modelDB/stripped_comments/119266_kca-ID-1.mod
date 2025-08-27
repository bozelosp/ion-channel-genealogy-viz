NEURON {
	SUFFIX kca
	USEION k READ ek WRITE ik
	USEION cal READ cali
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
	cali             (mM)           
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
	ik = gbar*po*(v + 85)    

}

DERIVATIVE states {
	rates(cali)
	m' = (m_inf - m) / tau_m
}


INITIAL {
	rates(cali)
	m = m_inf
}





PROCEDURE rates(cali(mM)) {  LOCAL a 

	a = cali/b
 	m_inf = a/(a+1)
	tau_m = taumin+ 1(ms)*1(mM)/(cali+b)
}