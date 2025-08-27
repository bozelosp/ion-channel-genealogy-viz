NEURON {
	SUFFIX KdrOlmKop
	USEION k READ ek WRITE ik
}
	
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {
    gkdr =   23 (mS/cm2)
    
}
    
ASSIGNED {
    ek      (mV)
    v       (mV)
    ik      (mA/cm2)
	ninf    (1)
	taon    (ms)
}

STATE { n }

INITIAL { 
    rates(v)
    n  = ninf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ik = (1e-3) * gkdr * n^4 * (v-ek)
}


DERIVATIVE states { 

    rates(v)
    n' = (ninf-n)/taon
}

PROCEDURE rates(v(mV)) { LOCAL an, bn

    an = fun3(v,   25,  -0.018,  -25)
    bn = fun3(v,   35,   0.0036,  12)
    
    ninf = an/(an+bn)
    taon = 1./(an+bn)
}

INCLUDE "custom_code/inc_files/135902_aux_fun.inc"