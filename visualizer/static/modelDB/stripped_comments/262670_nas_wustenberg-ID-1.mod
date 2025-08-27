INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
        SUFFIX nas
        USEION na READ ena WRITE ina
        RANGE gbar, ina, g
}

UNITS {
        (S) = (siemens)
        (mV) = (millivolt) 
        (mA) = (milliamp) 
}
 
PARAMETER { 
        gbar = 0.0      (mho/cm2)
}
 
ASSIGNED { 
	ena	(mV)
        v	(mV)
        ina	(mA/cm2)
        g	(S/cm2)
        minf
	hinf
        mtau	(ms)
        htau	(ms)
}
 
STATE {
    m
    h
}

BREAKPOINT { 
        SOLVE states METHOD cnexp 
        g = gbar * m * m * m * h
        ina = g * ( v - ena )
}
 
INITIAL { 
        settables(v)
	m = minf
        h  = hinf
} 

DERIVATIVE states { 
        settables(v) 
        h' = (hinf - h) / htau
	m' = (minf - m ) / mtau
}
















PROCEDURE settables(v (mV)) { 
UNITSOFF
        TABLE minf, hinf, mtau, htau FROM -120 TO 40 WITH 641
        minf  = 1.0 / (1 + exp((-30.1 - v) / 6.65))
        hinf  = 1.0 / (1 + exp((v + 51.4) / 5.9 ))
        mtau = (0.83 - 0.093) / (1 + exp((v + 20.3) / 6.45)) + 0.093
        htau = (12.24 - 1.9) / (1 + exp((v + 32.6) / 8.0)) + 1.9
UNITSON
}