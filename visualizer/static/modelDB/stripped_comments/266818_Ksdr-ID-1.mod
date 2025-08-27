NEURON {
    SUFFIX Ksdr 
    USEION k READ ek WRITE ik 
    
    RANGE g_Ks                              
    
    
    
    
    
    RANGE alpha_n                           
    
    RANGE ik                              
    
    RANGE beta_n                            
    
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
    
    g_Ks = 0.00000575 (mho/cm2) 
    v (mV)
    
    ek (mV)
    tscale = 0.001 (kHz)
}

ASSIGNED {
    
    alpha_n (1)                                
    
    ik  (mA/cm2)                                  
    
    beta_n (1)                                
    rate_nKs (/ms)
    
}

STATE {
    nKs  (1)
    
}

INITIAL {
    rates(v)
    rates(v) ? To ensure correct initialisation.
    
    nKs = 2.62753E-4
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    ik = ((  g_Ks   * (  nKs^2)) * (  v   -   ek  )) ? evaluable
}

DERIVATIVE states {
    rates(v)
    nKs' = rate_nKs 
    
}

PROCEDURE rates(v (mV)) {
    
    alpha_n = (4.81333E-6 (/mV) * ((  v   + 26.5 (mV)) / (1.0 - exp((- 0.128 (/mV) * (  v   + 26.5 (mV))))))) ? evaluable
    beta_n = (9.53333E-5 * exp(((- 0.038 (/mV)) * (  v   + 26.5 (mV))))) ? evaluable
    rate_nKs = tscale  * (((  alpha_n   * (1.0 -   nKs  )) - (  beta_n   *   nKs  ))) ? Note units of all quantities used here need to be consistent!
    
     
    
}