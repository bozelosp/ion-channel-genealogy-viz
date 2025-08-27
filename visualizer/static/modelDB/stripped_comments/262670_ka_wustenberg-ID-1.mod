INDEPENDENT { t FROM 0 TO 1 WITH 1 (ms) }

NEURON { 
        SUFFIX ka
        USEION k READ ek WRITE ik
        RANGE gbar, ik, g
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
	ek	(mV)
        v	(mV)
        ik	(mA/cm2)
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
        ik = g * ( v - ek )
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
        minf  = 1.0 / (1 + exp((-20.1 - v)/16.1))
        hinf  = 1.0 / ( 1 + exp( ( v + 74.7 ) / 7 ) )
        mtau = (1.65 - 0.35) / ((1 + exp(- (v + 70) / 4.0)) * (1 + exp((v + 20) / 12.0))) + 0.35
        htau = (90 - 2.5) / ((1 + exp(- (v + 60) / 25.0)) * (1 + exp((v + 62) / 16.0))) + 2.5
UNITSON
}