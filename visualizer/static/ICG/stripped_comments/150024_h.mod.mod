NEURON { 
	SUFFIX h
	NONSPECIFIC_CURRENT i
	RANGE gbar, m
	GLOBAL qdeltat, eh
} 
 
UNITS { 
	(mA) = (milliamp) 
	(mV) = (millivolt) 
} 
 
PARAMETER { 
    qdeltat = 1
    gbar = 1e-5 (siemens/cm2)
} 

ASSIGNED {
	v (mV)
	i (mA/cm2) 
	eh (mV)
	minf
	taum (ms)
} 
 
STATE {
	m
} 

INITIAL { 
    rate(v) 
	taum = 400 / qdeltat
    m = minf 
} 
 
BREAKPOINT { 
    SOLVE states METHOD cnexp 
	i = gbar * m * m * (v - eh)
} 
 
DERIVATIVE states { 
	rate(v) 
	m' = (minf - m) / taum 
} 

PROCEDURE rate(v(mV)) {
	TABLE minf FROM -150 TO 100 WITH 300 
	minf = 1 / (1 + exp((v + 80) / 5))
}