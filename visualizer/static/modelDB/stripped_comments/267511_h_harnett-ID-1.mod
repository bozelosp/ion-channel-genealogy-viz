UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    v                       (mV)
    celsius                 (degC)
    erev     = -30          (mV)
    gbar     = 0.0001       (mho/cm2)
    vhalf    = -100.6       (mV)
    k        = 6.4
    bA       = 9.63 
    bD       = 1.30 
    mA       = 0.0458 
    mD       = -0.0447 
    q10      = 2.2
    taumin	= 2.0	(ms)		
}

NEURON {
    SUFFIX h15
    NONSPECIFIC_CURRENT i
    RANGE gbar, minf, tau, g, m
    GLOBAL taumin, k, bA, bD, mA, mD, vhalf
}

STATE {
    m
}

ASSIGNED {
    i       (mA/cm2)
    minf
    tau
    g
}

INITIAL {
    rate(v)
    m = minf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g       = gbar*m
    i       = g*(v-erev)
}

DERIVATIVE states {     
    rate(v)
    m' = (minf - m) / tau
}

PROCEDURE rate(v (mV)) { 
    LOCAL qt
    qt = q10^((celsius-26.0)/10.0)

    if(v <= -92.0046199111992) {
      tau = exp(bA + mA * v) / qt 
    } else {
      tau = exp(bD + mD * v) / qt
    }
    if(tau < taumin) { tau = taumin }

    minf = 1.0/(1.0 + exp((v-vhalf)/k))
}