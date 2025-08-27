NEURON { 
	SUFFIX dcnCaHVA 
	USEION ca READ cai, cao WRITE ica
	RANGE perm, ica, m, cai
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
    perm = 7.5e-6 (cm/seconds)
} 

ASSIGNED {
	v (mV)
    cai (mM)
    cao (mM)     
	ica (mA/cm2) 
	minf (1)
	taum (ms) 
	celsius (degC)
	T (kelvin)
    A (1)
} 

STATE {
	m
} 

INITIAL { 
    T = 273.15 + celsius
    rate(v)
    m = minf 
} 
 
BREAKPOINT { 
    SOLVE states METHOD cnexp 
    A = getGHKexp(v)
    
    
    
    
    
    
    ica = perm * m*m*m * (4.47814e6 * v / T) * ((cai/1000) - (cao/1000) * A) / (1 - A)
} 
 
DERIVATIVE states { 
	rate(v) 
	m' = (minf - m)/taum 
} 

PROCEDURE rate(v(mV)) {
	TABLE minf, taum FROM -150 TO 100 WITH 300
	if(fabs(v + 34.5) < 1e-6) {
		minf = 1 / (1 + exp((v + 34.50001) / -9))
	} else {
	minf = 1 / (1 + exp((v + 34.5) / -9))
	}
    taum = 1 / ((31.746 * ((exp((v - 5) / -13.89) + 1) ^ -1)) + (3.97e-4 * (v + 8.9)) * ((exp((v + 8.9) / 5) - 1) ^ -1))
    taum = taum / qdeltat
} 

FUNCTION getGHKexp(v(mV)) {
    TABLE DEPEND T FROM -150 TO 100 WITH 300 
    getGHKexp = exp(-23.20764929 * v / T)
            
}