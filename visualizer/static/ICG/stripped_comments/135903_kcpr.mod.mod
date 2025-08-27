NEURON {

	SUFFIX kcpr
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gkc, ik
}
	
UNITS {    

    (mollar) = (1/liter)
	(mM)     = (millimollar)
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {

    gkc = 15    (mS/cm2)
    
}
    
ASSIGNED { 
    ek (mV)
    ik   (mA/cm2)    
    v    (mV)
    cai  (mM)
    cinf (1)
    tauc (ms)
}

STATE { c }

INITIAL { 
    
    rates(v)
    c  = cinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ik = (1e-3) * gkc * min(cai/250(mM),1) * c * (v-ek)
}


DERIVATIVE states { 

    rates(v)
    c' = (cinf-c)/tauc
}


PROCEDURE rates(v(mV)) { LOCAL a, b

    if (v<=-10) {
    
        a = 2(/ms) / 37.95 * ( exp( ( v + 50 ) / 11(mV) - ( v + 53.5 ) / 27(mV) ) )
        b = 2(/ms) * exp( ( - v - 53.5 ) / 27(mV) ) - a
    
    }else{
    
        a =  2(/ms) * exp( ( - v - 53.5 ) / 27(mV) )
        b = 0(/ms)
    
    }
    
    cinf = a/(a+b)
    tauc = 1.0/(a+b)    
}

INCLUDE "custom_code/inc_files/135903_aux_fun.inc"