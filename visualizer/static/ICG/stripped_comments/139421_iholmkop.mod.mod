NEURON {
	SUFFIX IhOlmKop
    NONSPECIFIC_CURRENT i
    RANGE gmax
    GLOBAL eh
}
	
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {
    gmax =  12    (mS/cm2)
    
}
    
ASSIGNED { 
    eh (mV)
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

FUNCTION qinf(v(mV))     { qinf = fun2(v, -84, 1, 10.2)*1(ms) }
FUNCTION qtau(v(mV))(ms) { qtau = 1(ms)/(exp((-14.59(mV)-0.086*v)/1(mV)) + exp((-1.87(mV)+0.0701*v)/1(mV))) }

INCLUDE "custom_code/inc_files/139421_aux_fun.inc"