INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX iKCa
    USEION k READ ek WRITE ik
    USEION ca READ cai
    RANGE gbar, m_inf, tau_m
    GLOBAL beta, cac
}


UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (molar) = (1/liter)
    (mM) = (millimolar)
}


PARAMETER {
    v                 (mV)
    celsius = 36      (degC)
    ek      = -80     (mV) 
    cai     = 50.0e-6 (mM)            
    gbar    = 0.002    (mho/cm2)
    beta    = 0.03    (1/ms)          
    
    cac     = 0.010   (mM)            
    taumin  = 0.1     (ms)            
}


STATE {m}        

ASSIGNED {       
    ik      (mA/cm2)
    m_inf
    tau_m   (ms)
    tadj
}
BREAKPOINT { 
    SOLVE states METHOD derivimplicit
    ik = gbar * m*m*m * (v - ek)    
}

DERIVATIVE states { 
    evaluate_fct(v,cai)
    m' = (m_inf - m) / tau_m
}

UNITSOFF
INITIAL {
    
    
    
    
    tadj = 3 ^ ((celsius-22.0)/10) 
    evaluate_fct(v,cai)
    m = m_inf
}

PROCEDURE evaluate_fct(v(mV),cai(mM)) {  LOCAL car
    car = (cai/cac)^2
    m_inf = car / ( 1 + car )      
    tau_m =  1 / beta / (1 + car) / tadj
    if(tau_m < taumin) { tau_m = taumin }   
}
UNITSON