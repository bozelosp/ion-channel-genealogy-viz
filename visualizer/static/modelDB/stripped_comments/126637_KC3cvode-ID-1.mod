UNITS {
	(molar) = (1/liter)
}

UNITS {
	(mV) =	(millivolt)
	(mA) =	(milliamp)
	(mM) =	(millimolar)
}


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX KC3cvode
	USEION ca READ cai
	USEION k WRITE ik
	RANGE gkbar,gk,ik
}


PARAMETER {
	celsius=37	(degC)
	v		(mV)
	gkbar=.08	(mho/cm2)	
	cai = .04e-3	(mM)
	ekcvode  = -85	(mV)
      mon = 1
	zon = 1	
}


ASSIGNED {
	ik		(mA/cm2)
	minf
	zinf
	gk
	tau 
	q10 
	alpha 
	beta 
	sum 
}

STATE {	m z }		

BREAKPOINT {
	SOLVE state METHOD cnexp

	ik = gkbar*m*z*z*(v - ekcvode)
}



INITIAL {
	
	m = minf
	z = zinf
}


DERIVATIVE state { 
	
	alpha = 400/(cai*1000)
	beta = 0.11/exp((v-35)/14.9)
	
	zinf = 1/(1+alpha)
      tau= 10
	z' = zon * (zinf-z)/tau      
	
	minf = 7.5/(7.5+beta)
      tau= 1/(7.5 + beta)
	m' = mon * (minf-m)/tau      

}