NEURON {
    SUFFIX Kti 
    USEION k READ ek WRITE ik 
    
    RANGE K1_leak_p                         
    RANGE K1_slope                          
    RANGE g                              
    RANGE K1_KIR_p                          
    
    
    
    
    
    RANGE ik                              
    
    RANGE K1_KIR_O                          
    
    RANGE K1_leak_O                         
    
}

UNITS {
    
    (nA) = (nanoamp)
    (uA) = (microamp)
    (mA) = (milliamp)
    (A) = (amp)
    (mV) = (millivolt)
    (mS) = (millisiemens)
    (uS) = (microsiemens)
    (kHz) = (kilohertz)
    (um) = (micrometer)
    (umol) = (micromole)
    (S) = (siemens)
    
}

PARAMETER {
    
    K1_leak_p = 0.5 (1)
    K1_slope = 11.16 (mV)
    g = 0.0002647 (mho/cm2) 
    K1_KIR_p = 0.5 (1)
    v (mV) 
    ek (mV) 
    tscale = 0.001 (kHz)
}

ASSIGNED {
    
    ik (mA/cm2)                                 
    
    K1_KIR_O   (mV)                            
    
    K1_leak_O     (mV)                         
    
}

INITIAL {
    rates(v)
    rates(v) ? To ensure correct initialisation.
    
}

BREAKPOINT {
    
    rates(v)
    ik = (  g   * (  K1_leak_O  +   K1_KIR_O  )) ? evaluable
    
    
    
}

PROCEDURE rates(v (mV)) {
    
    
    K1_KIR_O = ((  K1_KIR_p    * (  v   -   ek  )) / (1.0 + exp((  v   -   ek  ) /   K1_slope  ))) ? evaluable
    K1_leak_O = (  K1_leak_p   * (  v   -   ek  )) ? evaluable
    
     
    
}