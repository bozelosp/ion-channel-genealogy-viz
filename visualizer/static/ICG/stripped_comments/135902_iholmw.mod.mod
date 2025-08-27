NEURON {
	SUFFIX Iholmw
    NONSPECIFIC_CURRENT i
    GLOBAL eh
}
	
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {
    gh = 0.15 (mS/cm2)
    
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
	
    i = (1e-3) * gh * q * (v-eh)
}

DERIVATIVE states { q' = (qinf(v)-q)/qtau(v) }

FUNCTION qinf(v(mV))     { qinf = fun2(v, -80, 1, 10)*1(ms) }
FUNCTION qtau(v(mV))(ms) { qtau = 200(ms)/(exp((v+70(mV))/20(mV))+exp(-(v+70(mV))/20(mV))) + 5(ms) }

INCLUDE "custom_code/inc_files/135902_aux_fun.inc"