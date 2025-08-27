NEURON {
    SUFFIX NaK_pump 
    USEION na READ nai, nao WRITE ina 
    USEION k READ ki, ko WRITE ik 
    
    RANGE Km_Nai                            
    RANGE Km_Ko                             
    RANGE i_NaK_max                         
    
    
    
    
    
    
    
    
    
    RANGE f_NaK                             
    
    RANGE i_NaK                             
    
    RANGE ina
    RANGE ik
    RANGE sigma                             
    
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
    
    Km_Nai = 21000 (mM)
    Km_Ko = 1500 (mM)
    i_NaK_max = 0.88e-3 (mA/cm2)
    v (mV)
    nai (mM) 
    ki (mM)
    ko (mM)
    nao (mM) 
    tscale = 0.001 (kHz)
    celsius (degC)
}

ASSIGNED {
    
    T (degC)
    
    f_NaK (1)                                 
    
    i_NaK (mA/cm2)                                  
    
    ina (mA/cm2)
    
    ik (mA/cm2)
    
    sigma   (1)    
    
}


INITIAL {
    T = 273 + celsius
    
    rates(v)
    rates(v) ? To ensure correct initialisation.
    
}

BREAKPOINT {
    
    rates(v)
    i_NaK = (((  i_NaK_max   *   f_NaK  ) * (1.0 / (1.0 + ((  Km_Nai   /   nai  ) ^ 1.5)))) * (  ko   / (  ko   +   Km_Ko  ))) ? evaluable
    
    ina = 3 * i_NaK
    ik = -2 * i_NaK
}

PROCEDURE rates(v (mV)) {
    
    f_NaK = (1.0 / ((1.0 + (0.1245 * exp(((( - 0.1 ) *   v  ) * (  F   / (  R   *   T  )))))) + ((0.0365 *   sigma  ) * exp((( -   v  ) * (  F   / (  R   *   T  ))))))) ? evaluable

    sigma = ((1.0 / 7.0) * (exp((  nao   / 67300.0 (mM))) - 1.0)) ? evaluable
    
     
    
}