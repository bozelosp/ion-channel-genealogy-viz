NEURON {
    SUFFIX kca
    USEION k READ ek WRITE ik
    USEION ca READ cai
    RANGE gk, gbar, m_inf, tau_m, gmax
    GLOBAL beta, cac
}

UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (molar) = (1/liter)
    (mM) = (millimolar)
}

PARAMETER {
    v                (mV)
    celsius = 36     (degC)
    ek      = -80    (mV)
    cai     = 2.4e-5 (mM)           
    gbar    = 0.01   (mho/cm2)
    beta    = 0.03   (1/ms)          
    cac     = 0.025  (mM)            
    taumin  = 0.5    (ms)            
}

STATE {m}        

ASSIGNED {       
    ik      (mA/cm2)
    m_inf
    tau_m   (ms)
    tadj
    gk      (mho/cm2)
    gmax    (mho/cm2)
}
BREAKPOINT { 
    SOLVE states METHOD derivimplicit
    gk = gbar*m*m*m     
    ik = gk*(v - ek)    
    if (gk > gmax) {
        gmax = gk
    }
}

DERIVATIVE states { 
    evaluate_fct(v,cai)
    m' = (m_inf - m) / tau_m
}

INITIAL {
    
    
    
    
    tadj = 3 ^ ((celsius-22.0(degC))/10(degC)) 
    evaluate_fct(v,cai)
    m = m_inf
    gk = gbar*m*m*m
    gmax = gk
}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL car
    car = (cai/cac)^2
    m_inf = car / ( 1 + car )      
    tau_m =  1 / beta / (1 + car) / tadj
    if(tau_m < taumin) { tau_m = taumin }   
}