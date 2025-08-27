NEURON {
	SUFFIX kcain
	USEION k READ ko, ki WRITE ik
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
	ki 		(mM)
	ko		(mM)
	m_inf
	tau_m           (ms)



}

PARAMETER {
	gbar    = 10   (mho/cm2)

	taumin  = 0	(ms)  
	b 	= 0.008 (/ms)  
	

}


STATE {
	m   
}

BREAKPOINT { 
	SOLVE states METHOD cnexp
	ek = 25 * log(ko/ki)
	po = m*m
	ik = gbar*po*(v - ek)    
}

DERIVATIVE states {
	rates(cai)


	m' = (m_inf - m) / tau_m 

	
} 


INITIAL {
	rates(cai)
	m = 0



}


PROCEDURE rates(cai(mM)) { 
	LOCAL a



	

	a = cai/b
	m_inf = a/(a+1)

	tau_m = taumin+ 1(ms)*1(mM)*b/(cai+b)



}