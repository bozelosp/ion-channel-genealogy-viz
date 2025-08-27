UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
    (S) = (siemens)
    (um) = (micrometer)
    (molar) = (1/liter)
    (mM) = (millimolar)
    (l) = (liter)
}


    
NEURON {
      

    SUFFIX Leak
    
    RANGE e
    NONSPECIFIC_CURRENT il
    
    RANGE gmax, gion,il
    
}

PARAMETER { 
      

    gmax = 0.0003 (S/cm2) 
    
    e = -80 (mV) 
    
}



ASSIGNED {
      

    v (mV)
        
    il (mA/cm2)
        
}

BREAKPOINT { 
    il = gmax*(v - e) 
        

}