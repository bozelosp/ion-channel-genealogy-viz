NEURON {
	SUFFIX IhPyrKop
    NONSPECIFIC_CURRENT i
	RANGE v50, gmax
}
	
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {
    gmax =   0.0  (mS/cm2)
    eh   = -30.0  (mV)
    v50  =   0.0 (mV)
}
    
ASSIGNED { 

    v (mV)
    i (mA/cm2)
}

STATE { q }

INITIAL { q  = qinf(v) }

BREAKPOINT {

	SOLVE states METHOD cnexp
	
    i = (1e-3) * gmax * q * (v-eh)
}

DERIVATIVE states { q' = (qinf(v)-q)/qtau(v) }

FUNCTION qinf(v(mV))     { qinf = fun2(v, v50, 1.0, 10.5)*1.0(ms) }
FUNCTION qtau(v(mV))(ms) { qtau = 1.0(ms)/(exp((-14.59(mV)-0.086*v)/1.0(mV)) + exp((-1.87(mV)+0.0701*v)/1.0(mV))) }





INCLUDE "aux_fun.inc"