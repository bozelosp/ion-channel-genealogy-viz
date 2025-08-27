UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
     (mM) = (milli/liter)

}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {
	dt 	   		(ms)
	v 	   		(mV)
        ehd=-47 		(mV) 				       
	gbar=0 (pS/um2)	
	gamma_ih	
	seed		
	vshift = 0
}


NEURON {
	SUFFIX ih
	NONSPECIFIC_CURRENT Iqq
	RANGE Iqq,gbar,vshift,ehd, qtau, qinf, gq
}

STATE {
	qq
}

ASSIGNED {
	Iqq (mA/cm2)
	qtau (ms)
	qinf
	gq	(pS/um2)
	
}

INITIAL {
	qq=alpha(v-vshift)/(beta(v-vshift)+alpha(v-vshift))

	qtau = 1./(alpha(v-vshift) + beta(v-vshift))
	qinf = alpha(v)/(alpha(v-vshift) + beta(v-vshift))
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	
	qtau = 1./(alpha(v-vshift) + beta(v-vshift))
	qinf = alpha(v-vshift)/(alpha(v-vshift) + beta(v-vshift))
	
	gq = gbar*qq
	Iqq = (1e-4)*gq*(v-ehd)
	
}

FUNCTION alpha(v(mV)) {

	alpha = 0.001*6.43*(v+154.9)/(exp((v+154.9)/11.9)-1)
	
        
        

}

FUNCTION beta(v(mV)) {
	beta = 0.001*193*exp(v/33.1)			
}

DERIVATIVE state {     
	qq' = (1-qq)*alpha(v-vshift) - qq*beta(v-vshift)
}