INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
THREADSAFE
	SUFFIX iahp
	USEION k READ ek WRITE ik VALENCE 1
	USEION Ca READ Cai VALENCE 2
      RANGE gkbar, g, minf, taum
	GLOBAL beta, cac, m_inf, tau_m, x
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}


PARAMETER {
	v (mV)
	ek = -90 (mV)
	celsius = 36 (degC)
	Cai 	= 5e-5 (mM)			
	gkbar	= 1.3e-4	(mho/cm2)	
	beta	= 0.02	(1/ms)	
	cac	= 4.3478e-4(mM)		
	taumin	= 1	(ms)		
      x       = 2				
}




STATE {
	m
}


ASSIGNED {
	ik 	(mA/cm2)
	g       (mho/cm2)
	m_inf
	tau_m	(ms)
	minf
      taum
	tadj
}


BREAKPOINT { 
	SOLVE states METHOD cnexp
        minf = m_inf
        taum = tau_m
	  g = gkbar*m*m
	  ik = g * (v - ek)
}

DERIVATIVE states { 
	evaluate_fct(v,Cai)
	m' = (m_inf - m) / tau_m
}


UNITSOFF
INITIAL {



	VERBATIM
	Cai = _ion_Cai;
	ENDVERBATIM

	tadj = 3 ^ ((celsius-22.0)/10)
	evaluate_fct(v,Cai)
	m = m_inf
      minf = m_inf
      taum = tau_m
}

PROCEDURE evaluate_fct(v(mV),Cai(mM)) {  LOCAL car, tcar
	car = (Cai/cac)^x
	m_inf = car / ( 1 + car )
	tau_m = 1 / beta / (1 + car) / tadj
      if(tau_m < taumin) { tau_m = taumin } 	
}

UNITSON