NEURON {
    SUFFIX Kto_f 
    USEION k READ ki, ko WRITE ik 

    RANGE g_Kto_f                           
    
    
    
    
    RANGE ko                                
    RANGE ki                                
    
    
    
    RANGE beta_a                            
    
    RANGE alpha_a                           
    
    RANGE fast_transient_outward_potassium_current_alpha_i 
    
    RANGE ik                           
    
    
    
    RANGE fast_transient_outward_potassium_current_beta_i 
    
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
    F       = (faraday) (kilocoulombs)
    
    R       = (k-mole) (joule/degC)
    
}

PARAMETER {
    
    g_Kto_f = 0.0004067 (mho/cm2)
    
    v (mV)
    ko (mM)
    ki (mM)
    celsius (degC)
    
    
    
}

ASSIGNED {
    
    T (degC)
    
    beta_a (/ms)                                
    
    alpha_a (/ms)                               
    
    fast_transient_outward_potassium_current_alpha_i (/ms) 
    
    ik (mA/cm2)                                
    
    E_K (mV)                                   
    
    fast_transient_outward_potassium_current_beta_i (/ms) 
    rate_ato_f (/ms)
    rate_ito_f (/ms)
    
}

STATE {
    ato_f  (1)
    ito_f  (1)
    
}

INITIAL {
    
    T = 273 + celsius

    rates(v)
    rates(v) ? To ensure correct initialisation.
    
    ato_f = 0.00265563
    
    ito_f = 0.999977
    
}

BREAKPOINT {
    
    SOLVE states METHOD cnexp
    ik = (((  g_Kto_f   * (  ato_f   ^ 3)) *   ito_f  ) * (  v   -   E_K  )) ? evaluable
    
    
}

DERIVATIVE states {
    rates(v)
    ato_f' = (((  alpha_a   * (1.0 -   ato_f  )) - (  beta_a   *   ato_f  ))) 
    ito_f' = (((  fast_transient_outward_potassium_current_alpha_i   * (1.0 -   ito_f  )) - (  fast_transient_outward_potassium_current_beta_i   *   ito_f  ))) 
    
}

PROCEDURE rates(v (mV)) {
    
    beta_a = (0.3956 (/ms) * exp(((- 0.06237(/mV)) * (  v   + 30.0(mV))))) ? evaluable
    alpha_a = (0.18064 (/ms) * exp((0.03577 (/mV)* (  v   + 30.0  (mV))))) ? evaluable
    fast_transient_outward_potassium_current_alpha_i = (1.52E-4 (/ms)* (exp((( - (  v   + 13.5 (mV))) / 7.0 (mV))) / ((0.0067083 * exp(((- (  v   + 33.5 (mV))) / 7.0 (mV)))) + 1.0))) ? evaluable
    E_K = ((  R   * (  T   /   F  )) * log((  ko   /   ki  ))) ? evaluable
    fast_transient_outward_potassium_current_beta_i = (9.5E-4 (/ms) * (exp(((  v   + 33.5 (mV)) / 7.0(mV))) / ((0.051335 * exp(((  v   + 33.5 (mV)) / 7.0 (mV)))) + 1.0))) ? evaluable
    
    
    
     
    
}