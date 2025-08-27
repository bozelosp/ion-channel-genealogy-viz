NEURON {
	SUFFIX KaOlmKop
	USEION k READ ek WRITE ik
	GLOBAL ek
}
	
UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(mS) = (millisiemens)
}

PARAMETER {
    gka =   16 (mS/cm2)
    ek  =  -90 (mV)
}
    
ASSIGNED {
    v       (mV)
    ik      (mA/cm2)
}

STATE { a b }

INITIAL { 
    a  = ainf(v)
    b  = binf(v)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = (1e-3) * gka * a * b * (v-ek)
}


DERIVATIVE states { 
	a' = (ainf(v)-a)/atau(v) 
	b' = (binf(v)-b)/btau(v) 
}

FUNCTION ainf(v(mV))     { ainf = fun2(v, -14, 1, -16.6)*1(ms) }
FUNCTION atau(v(mV))(ms) { atau = 5(ms) }

FUNCTION binf(v(mV))     { binf = fun2(v, -71, 1, 7.3)*1(ms) }
FUNCTION btau(v(mV))(ms) { btau = 1(ms)/(0.000009*exp(-(v-26)/18.5(mV)) + 0.014/(0.2+exp(-(v+70)/11(mV)))) }

INCLUDE "custom_code/inc_files/182134_aux_fun.inc"