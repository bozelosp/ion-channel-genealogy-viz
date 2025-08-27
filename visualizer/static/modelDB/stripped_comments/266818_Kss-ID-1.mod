NEURON {
    SUFFIX Kss 
    USEION k READ ek, ki, ko WRITE ik 
    RANGE g                             
    
    
    
    RANGE ass                               
    
    
    RANGE ik                             
    
    RANGE tau_Kss                           
    
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
    
    g = 0.000088 (mho/cm2)
    v (mV) 
    
    
    
}

ASSIGNED {
    
    ass (1) 
    ek (mV)                                   
    ik (mA/cm2) 
    ko (mM)
    ki (mM)
    tau_Kss (ms)                                
    rate_aKss (/ms)
    rate_iKss (/ms)
    
    iKss (1)
}

STATE {
    aKss (1)
    
    
}

INITIAL {
    rates(v)
    rates(v) ? To ensure correct initialisation.
    
    aKss = 4.17069E-4
    
    iKss = 1.0
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    ik = (((  g   *   aKss  ) *   iKss  ) * (  v   -   ek  )) ? evaluable
    
    
}

DERIVATIVE states {
    rates(v)
    aKss' = (((  ass   -   aKss  ) /   tau_Kss  ))
    
    
}

PROCEDURE rates(v (mV)) {
    
    ass = 1 / (1 + exp(-(v + 22.5 (mV) ) / 7.7 (mV) ))
    tau_Kss = ((39.3 (ms) * exp((( - 0.0862 (/mV)) *   v  ))) + 13.17 (ms)) ? evaluable
    
    
    
     
    
}