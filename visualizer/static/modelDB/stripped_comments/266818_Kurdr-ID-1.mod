NEURON {
    SUFFIX Kurdr 
    USEION k READ ek WRITE ik 
    RANGE g_Kur                             
    
    RANGE iss                               
    
    
    RANGE ass                               
    
    
    RANGE tau_iur                           
    
    RANGE tau_aur                           
    
    RANGE ik                             
    
}

UNITS {
    
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (molar) = (1/liter)
    (kHz) = (kilohertz)
    (mM) = (millimolar)
    (um) = (micrometer)
    (umol) = (micromole)
    (S) = (siemens)
    
}

PARAMETER {
    
    g_Kur = 0.000112 (S/cm2)
    
    v (mV)
    ek (mV)
    
}

ASSIGNED {
    
    tau_iur (ms)                               
    
    tau_aur (ms)                               
    ass (1) 
    iss (1) 
    
    ik (mA/cm2)                                 
    rate_aur (/ms)
    rate_iur (/ms)
    
}

STATE {
    aur  (1)
    iur  (1)
    
}

INITIAL {
    rates(v)
    rates(v) ? To ensure correct initialisation.
    
    aur = 4.17069E-4
    
    iur = 0.998543
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    ik = (((  g_Kur   *   aur  ) *   iur  ) * (  v   -   ek  )) ? evaluable
}

DERIVATIVE states {
    rates(v)
    aur' =  (((  ass   -   aur  ) /   tau_aur  ))
    iur' = (((  iss   -   iur  ) /   tau_iur  ))
    
}

PROCEDURE rates(v (mV)) {
    
    ass = 1 / (1 + exp(-(v + 22.5 (mV) ) / 7.7 (mV) ))
    iss = 1 / (1 + exp((v + 45.2 (mV)) / 5.7 (mV)))
    tau_iur = (1200.0 (ms)- (170.0 (ms)/ (1.0 + exp(((  v   + 45.2 (mV)) / 5.7(mV)))))) ? evaluable
    tau_aur = ((0.493 (ms)* exp(((- 0.0629 (/mV)) *   v  ))) + 2.058(ms)) ? evaluable
    
    
}