NEURON {
    SUFFIX ih
    NONSPECIFIC_CURRENT i
    RANGE i, gslow, gfast, gslowbar, gfastbar
    GLOBAL ehcn, taufn, taufdo, taufdd, taufro, taufrd
    GLOBAL tausn, tausdo, tausdd, tausro, tausrd
    GLOBAL mifo, mifd, mife, miso, misd, mise
}

UNITS {
    (mV) = (millivolt)
    (S) = (siemens)
    (mA) = (milliamp)
}

PARAMETER {
    gfastbar = 9.8e-5    (S/cm2)
    gslowbar = 5.3e-5    (S/cm2)
    ehcn    = -20        (mV)
    taufn   = 0.51       (ms)    
    taufdo  = 1.7        (mV)
    taufdd  = 10         (mV)
    taufro    = 340      (mV)
    taufrd    = 52       (mV)
    tausn   = 5.6        (ms)    
    tausdo  = 17         (mV)
    tausdd  = 14         (mV)
    tausro    = 260      (mV)
    tausrd    = 43       (mV)
    mifo    = 74.2       (mV)    
    mifd    = 9.78       (mV)
    mife    = 1.36
    miso    = 2.83       (mV)    
    misd    = 15.9       (mV)
    mise    = 58.5
}

ASSIGNED {
    v        (mV)
    gslow    (S/cm2)
    gfast    (S/cm2)
    i        (mA/cm2)
    alphaf   (/ms)        
    betaf    (/ms)        
    alphas   (/ms)        
    betas    (/ms)        
}

INITIAL {
    
    settables(v)
    mf = alphaf/(alphaf+betaf)
    ms = alphas/(alphas+betas)
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    gfast = gfastbar*mf
    gslow = gslowbar*ms
    i = (gfast+gslow)*(v-ehcn)
}

STATE {
    mf ms
}

DERIVATIVE states {  
    settables(v)
    mf' = alphaf*(1-mf) - betaf*mf
    ms' = alphas*(1-ms) - betas*ms
}

PROCEDURE settables(v (mV)) { 
    LOCAL mif, mis, tauf, taus 
    TABLE alphaf, betaf, alphas, betas FROM -100 TO 100 WITH 200

    tauf = taufn/( exp( (v-taufdo)/taufdd ) + exp( -(v+taufro)/taufrd ) )
    taus = tausn/( exp( (v-tausdo)/tausdd ) + exp( -(v+tausro)/tausrd ) )
    mif = 1/pow( 1 + exp( (v+mifo)/mifd ), mife )
    mis = 1/pow( 1 + exp( (v+miso)/misd ), mise )

    alphaf = mif/tauf
    alphas = mis/taus
    betaf = (1-mif)/tauf
    betas = (1-mis)/taus
}