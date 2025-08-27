NEURON {
    SUFFIX NAbg 
    USEION na READ ena WRITE ina
    
    RANGE g                             
    
    
    
    
    RANGE ina                             
    
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
    
    g = 0.0000029 (mho/cm2) 
    v (mV)
    ena (mV) 
    tscale = 0.001 (kHz)
}

ASSIGNED {
    
    ina (mA/cm2)                                  
    
}


INITIAL {
    rates(v)
    rates(v) ? To ensure correct initialisation.
    
}

BREAKPOINT {
    
    rates(v)
    
}

PROCEDURE rates(v (mV)) {
    
    ina = (  g   * (  v   -   ena  )) ? evaluable
    
     
    
}