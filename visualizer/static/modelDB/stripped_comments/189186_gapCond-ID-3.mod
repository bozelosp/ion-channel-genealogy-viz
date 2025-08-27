?  This is a NEURON mod file generated from a ChannelML file

?  Unit system of original ChannelML file



? Creating synaptic mechanism for an electrical synapse
    






UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

    
NEURON {
    POINT_PROCESS gapCond
    NONSPECIFIC_CURRENT i
    RANGE g, i
    RANGE weight
    
    POINTER vgap  
        

}

PARAMETER {
    v (millivolt)
    vgap (millivolt)
    g = 1e-6 (microsiemens)
    weight = 1

}


ASSIGNED {
    i (nanoamp)
}

BREAKPOINT {
    i = weight * g * (v - vgap)
}