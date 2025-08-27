NEURON {
	SUFFIX NafOlmKop
	USEION na READ ena WRITE ina
}
	
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {
    gna  = 30 (mS/cm2)
    
}
    
ASSIGNED {
    ena     (mV)
    v       (mV)
    ina     (mA/cm2)
    minf    (1)
    hinf    (1)
    taom    (ms)
    taoh    (ms)
}

STATE { m h }

INITIAL {
    rates(v)
    m    = minf
    h    = hinf
}

BREAKPOINT {

	SOLVE states METHOD cnexp
	
	ina = (1e-3) * gna * m^3 * h * (v-ena)
}

DERIVATIVE states { 
    rates(v)
    m' = (minf-m)/taom
    h' = (hinf-h)/taoh
}

PROCEDURE rates(v(mV)) { LOCAL am, bm, ah, bh
    
    am   = fun3(v,  -38, -0.1,    -10)
    bm   = fun1(v,  -65,  4,      -18)
    minf = am/(am+bm)
    taom = 1./(am+bm)
 
    ah   = fun1(v,  -63,    0.07,  -20)
    bh   = fun2(v,  -33,    1,     -10)
    hinf = ah/(ah+bh)
    taoh = 1./(ah+bh)
}

INCLUDE "custom_code/inc_files/135903_aux_fun.inc"