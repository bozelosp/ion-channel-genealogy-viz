NEURON { 
	SUFFIX CaLVA 
	
	USEION ca READ cai, cao WRITE ica
        RANGE perm, ica, m, h, cai
	GLOBAL qdeltat
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
} 
 
PARAMETER { 
    qdeltat = 1
    perm = 1 (cm/seconds)
} 

ASSIGNED {
    v (mV)
    cai (mM)
    cao (mM)     
	ica (mA/cm2) 
	minf
	hinf
	taum (ms) 
	tauh (ms) 
	celsius (degC)
	T (kelvin)
    A (1)
} 
 
STATE {
	m
    h
} 

INITIAL { 
    T = 273.15 + celsius
    rate(v)
    m = minf 
	h = hinf
} 
 
BREAKPOINT { 
    SOLVE states METHOD cnexp 
    A = getGHKexp(v)
    
    
    
    
    
    
    ica = perm * m*m * h * (4.47814e6 * v / T) * ((cai/1000) - (cao/1000) * A) / (1 - A)
} 
 
DERIVATIVE states { 
	rate(v) 
	m' = (minf - m)/taum 
	h' = (hinf - h)/tauh 
} 

PROCEDURE rate(v(mV)) {
	TABLE minf, taum, hinf, tauh  FROM -150 TO 100 WITH 300 
	minf = 1 / (1 + exp((v + 56) / -6.2))
	taum = 0.333 / (exp((v + 131) / -16.7) + exp((v + 15.8) / 18.2)) + 0.204
    taum = taum / qdeltat
	hinf = 1 / (1 + exp((v + 80) / 4))
    if (v < -81) {
        tauh = 0.333 * exp((v + 466) / 66)
    } else {
        tauh = 0.333 * exp((v + 21) / -10.5) + 9.32
    }
    tauh = tauh / qdeltat
}

FUNCTION getGHKexp(v(mV)) {
    TABLE DEPEND T FROM -150 TO 100 WITH 300 
    getGHKexp = exp(-23.20764929 * v / T)
            
}