NEURON {
    SUFFIX Ksto 
    USEION k READ ek WRITE ik 
    
    RANGE g_Kto_s                           
    
    
    
    
    
    RANGE ass                               
    
    RANGE ik                          
    
    RANGE iss                               
    
    RANGE tau_ti_s                          
    
    RANGE tau_ta_s                          
    
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
    
    g_Kto_s = 0 (mho/cm2)
    v (mV) 
    ek (mV) 
    tscale = 0.001 (kHz)
}

ASSIGNED {
    
    ass (1)                                  
    
    ik (mA/cm2)                                
    
    iss  (1)                                  
    
    tau_ti_s (ms)                              
    
    tau_ta_s (ms)                              
    rate_ato_s (/ms)
    rate_ito_s (/ms)
    
}

STATE {
    ato_s  (1)
    ito_s  (1)
    
}

INITIAL {
    rates(v)
    rates(v) ? To ensure correct initialisation.
    
    ato_s = 4.17069E-4
    
    ito_s = 0.998543
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    ik = (((  g_Kto_s   *   ato_s  ) *   ito_s  ) * (  v   -   ek  )) ? evaluable    
    
}

DERIVATIVE states {
    rates(v)
    ato_s' = rate_ato_s 
    ito_s' = rate_ito_s 
    
}

PROCEDURE rates(v (mV)) {
    
    ass = (1.0 / (1.0 + exp(((0.0 - (  v   + 22.5 (mV))) / 7.7 (mV))))) ? evaluable

    iss = (1.0 / (1.0 + exp(((  v   + 45.2(mV)) / 5.7(mV))))) ? evaluable
    tau_ti_s = (270.0 (ms)+ (1050.0 (ms) / (1.0 + exp(((  v   + 45.2(mV)) / 5.7(mV)))))) ? evaluable
    tau_ta_s = ((0.493 (ms) * exp((( - 0.0629(/mV)) *   v  ))) + 2.058(ms)) ? evaluable
    rate_ito_s = (((  iss   -   ito_s  ) /   tau_ti_s  )) ? Note units of all quantities used here need to be consistent!
    rate_ato_s = (((  ass   -   ato_s  ) /   tau_ta_s  )) ? Note units of all quantities used here need to be consistent!
    
     
    
}