NEURON {
	SUFFIX ICaolmw
	USEION ca WRITE ica
}
	
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {
    gca = 1    (mS/cm2)
    eca = 120  (mV)
}
    
ASSIGNED { 
    ica (mA/cm2)    
    v   (mV)
}

BREAKPOINT { ica = (1e-3) * gca * mcainf(v)^2 * (v-eca) }

FUNCTION mcainf(v(mV)) { mcainf = fun2(v, -20, 1,  -9)*1(ms) }

INCLUDE "aux_fun.inc"