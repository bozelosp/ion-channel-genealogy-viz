INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX iCaAN
	USEION can READ ecan WRITE ican VALENCE 1
	USEION ca READ cai
        RANGE gbar, m_inf, tau_m
	RANGE ican
	GLOBAL beta, cac, taumin
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


PARAMETER {
	v		  (mV)
	celsius		  (degC)
        dt                (ms)
	ecan	= -20	  (mV)		
	cai		  (mM)
	gbar	= 0.00025 (mho/cm2)
	beta	= 2.0e-3  (1/ms)	

	tau_factor = 40         

	
	cac		= 5e-4	  (mM)		

	taumin	= 0.1	  (ms)		
	
	
	
	
	
	
	
}


STATE {
	m
}

ASSIGNED {
	ican	(mA/cm2)
	m_inf
	tau_m	(ms)
	tadj
}

BREAKPOINT { 
	SOLVE states 
	ican = gbar * m*m * (v - ecan)
}






  
PROCEDURE states() {
        evaluate_fct(v,cai)
	
        m = m + ( 1-exp(-dt/tau_m) )*(m_inf-m)
	

}

UNITSOFF
INITIAL {




	tadj = 3.0 ^ ((celsius-22.0)/10)

	evaluate_fct(v,cai)
	m = m_inf
}


PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL alpha2

	alpha2 = beta * (cai/cac)^2
	
	tau_m = tau_factor / (alpha2 + beta) / tadj		
	
	m_inf = alpha2 / (alpha2 + beta)							

	if(tau_m < taumin) { tau_m = taumin }					

}
UNITSON