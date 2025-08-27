NEURON {
	SUFFIX KahpM95
	USEION k READ ek WRITE ik
        USEION ca READ cai
        RANGE  gbar,gkahp,ik
        GLOBAL inf,tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	celsius = 6.3	(degC)
	gbar	= .003 	(mho/cm2)
        n	= 4
        cai	= 50.e-6 (mM)
        a0	= 1.3e13 (/ms-mM-mM-mM-mM)	
        b0	= .5e-2  (/ms)			
        v       	 (mV)
        ek      	 (mV)
}

STATE {	w }

ASSIGNED {
	ik 		(mA/cm2)
        gkahp  		(mho/cm2)
        inf
        tau
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gkahp = gbar*w
	ik = gkahp*(v-ek)
}

INITIAL {
	rate(cai)
	w=inf
}

FUNCTION alp(cai (mM)) {
  alp = a0*cai^n
}

DERIVATIVE state {     
        rate(cai)
        w' = (inf - w)/tau
}

PROCEDURE rate(cai (mM)) { 
        LOCAL a
        a = alp(cai)
        tau = 1/(a + b0)
        inf = a*tau
}