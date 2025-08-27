UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
}

PARAMETER {
    v                       (mV)
    celsius  = 34.0    (degC)
    erev     = -37.0               (mV)
    gbar     = 0.0001       (mho/cm2)
    vhalfl   = -78.474      (mV)    
    kl       = -6                   
    vhalft   = -66.139      (mV)    
    a0t      = 0.009        (/ms)   
    zetat    = 20           (1)     
    gmt      = 0.01         (1)     
    q10      = 4.5
    qtl      = 1
    taumin	= 2.0	(ms)		
}

NEURON {
    SUFFIX h
    NONSPECIFIC_CURRENT i
    RANGE gbar, vhalfl
    RANGE linf, taul, g
    GLOBAL taumin
}

STATE {
    l
}

ASSIGNED {
    i       (mA/cm2)
    linf
    taul
    g
}

INITIAL {
    rate(v)
    l       = linf
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g       = gbar*l
    i       = g*(v-erev)
}

FUNCTION alpt(v(mV)) {
    alpt    = exp(0.0378*zetat*(v-vhalft))
}

FUNCTION bett(v(mV)) {
    bett    = exp(0.0378*zetat*gmt*(v-vhalft))
}

DERIVATIVE states {     
    rate(v)
    l'      = (linf - l)/taul
}

PROCEDURE rate(v (mV)) { 
    LOCAL a,qt
    qt      = q10^((celsius-33)/10)
    a       = alpt(v)
    linf    = 1/(1 + exp(-(v-vhalfl)/kl))
    taul    = bett(v)/(qtl*qt*a0t*(1+a)) + 1e-8 
    if(taul < taumin) { taul = taumin } 	
}