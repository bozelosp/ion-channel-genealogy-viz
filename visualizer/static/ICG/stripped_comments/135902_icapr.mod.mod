NEURON {

	SUFFIX icapr
	USEION ca READ eca WRITE ica
	RANGE gca, ica
}
	
UNITS {

	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {

    gca = 10  (mS/cm2)
    
}
    
ASSIGNED { 
    eca  (mV)
    ica  (mA/cm2)    
    v    (mV)
    sinf (1)
    taus (ms)
}

STATE { s }

INITIAL { 

    rates(v)
    s  = sinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ica = (1e-3) * gca * s^2 * (v-eca)
}


DERIVATIVE states { 

    rates(v)
    s' = (sinf-s)/taus
}

PROCEDURE rates(v(mV)) { LOCAL a,b

    a = fun2(v,      5,    1.6,     -1/0.072)
    b = fun3(v,   -8.9,   0.02,     5)

    sinf = a/(a+b)
    taus = 1.0/(a+b)
}

INCLUDE "custom_code/inc_files/135902_aux_fun.inc"