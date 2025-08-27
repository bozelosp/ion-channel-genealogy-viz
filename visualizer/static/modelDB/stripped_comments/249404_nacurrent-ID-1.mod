UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
}
 
NEURON {
    SUFFIX nacurrent
    USEION na READ ena WRITE ina
    RANGE gnabar, gna
    GLOBAL minf, hinf, mtau, htau
    THREADSAFE
}
 
PARAMETER {
    gnabar = .0 (S/cm2)	<0, 1e9>
}
 
STATE {
    m (1)
    h (1)
}
 
ASSIGNED {
    v (mV)
    celsius (degC)
    ena (mV)
    gna (S/cm2)
    ina (mA/cm2)
    minf (1)
    hinf (1)
    mtau (ms)
    htau (ms)
}
 
BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gnabar * m * m * m * h
    ina = gna * (v - ena)
}
 
INITIAL {
    rates(v)
    m = minf
    h = hinf
}

DERIVATIVE states {  
    rates(v)
    m' = (minf - m) / mtau
    h' = (hinf - h) / htau
}
 
PROCEDURE rates(v(mV)) {
    
    LOCAL  alpha, beta, sum, q10
    TABLE minf, mtau, hinf, htau DEPEND celsius FROM -100 TO 100 WITH 200

UNITSOFF
    q10 = 3 ^ ((celsius - 23) / 10)
    
    alpha = -.182 * vtrap(-(v + 40), 6)
    beta =  -.124 * vtrap((v + 40), 6)
    sum = alpha + beta
    mtau = 0.25 / (q10 * sum)
    minf = alpha / sum
    
    alpha = -0.015 * vtrap((v + 66), 6)
    beta = -0.015 * vtrap(-(v + 66), 6)
    sum = alpha + beta
    htau = 1 / (q10 * sum)
    hinf = alpha / sum
}
 
FUNCTION vtrap(x, y) {
    
    if (fabs(x / y) < 1e-6) {
        vtrap = -y * (1 - x / y / 2)
    } else {
        vtrap = x / (1 - exp(x / y))
    }
}
 
UNITSON