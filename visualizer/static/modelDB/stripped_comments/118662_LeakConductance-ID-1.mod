?  This is a NEURON mod file generated from a v1.7.2 ChannelML file

?  Unit system of original ChannelML file








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
      

    SUFFIX LeakConductance
    ? A non specific current is present
    RANGE e
    NONSPECIFIC_CURRENT i
    
    RANGE gmax, gion
    
}

PARAMETER { 
      

    gmax = 0.0003 (S/cm2) 
    
    e = -60 (mV)
    
}



ASSIGNED {
      

    v (mV)
        
    i (mA/cm2)
        
}

BREAKPOINT { 
    i = gmax*(v - e) 
        

}