NEURON {
	SUFFIX icnew
	USEION k READ ek WRITE ik
	USEION ca READ cai
	GLOBAL inf
	RANGE ica,gkbar,ek,tau_diff,taum,fact, ik, ca_beta
}
UNITSON
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
PARAMETER {
	v (mV)
	dt (ms)
	gkbar=.045 (mho/cm2)
	ek = -95 (mV)
	
        tau_diff = 0.05 (ms)
	taum = 2.0 (1/ms)
	fact = 1e-3

}

CONSTANT {
	ca_beta  = 20.0 (1/ms)   
	ca_alpha = 100.0 (mM/ms/mA)
}

STATE { m}
ASSIGNED {
	ik (mA/cm2)
	ica (mA/cm2)
	cai (mM)
	inf[1]
	tau[1]
}

INITIAL {
	
  	
	
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gkbar*m*m*(v - ek)
}

DERIVATIVE states {	
	mhn(cai)
	m' = (inf[0] - m)/tau[0]
	
}

FUNCTION varss(ca) {
	varss = ca / (ca + 0.040) 
	
}

FUNCTION vartau() {
	vartau = taum  
}	


PROCEDURE mhn(ca) { 



	tau[0] = vartau()
	inf[0] = varss(ca)
}