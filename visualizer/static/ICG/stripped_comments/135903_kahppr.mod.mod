NEURON {

	SUFFIX kahppr
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gkahp, ik, qinf, tauq
}
	
UNITS {

    (mollar) = (1/liter)
	(mM)     = (millimollar)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {

    gkahp = 0.8 (mS/cm2)
    
}
    
ASSIGNED {
    ek (mV)
    v    (mV)
    ik   (mA/cm2)
    cai  (mM)
    qinf (1)
    tauq (ms)
}

STATE { q }

INITIAL { 
    
    rates(v)
    q  = qinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ik = (1e-3) * gkahp * q * (v-ek)
}


DERIVATIVE states { 

    rates(v)
    q' = (qinf-q)/tauq
        
}

PROCEDURE rates(v(mV)) { LOCAL a,b
    
    a = 0.01(/ms) * min(cai/500(mM),1)
    b = 1(/ms)/1000
    
    qinf = a/(a+b)
    tauq = 1.0/(a+b)
}

INCLUDE "custom_code/inc_files/135903_aux_fun.inc"