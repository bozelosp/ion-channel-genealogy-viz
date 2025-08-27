UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}



NEURON {
    THREADSAFE
    
    SUFFIX h

    NONSPECIFIC_CURRENT i

    
    
    
    RANGE gmax, tau, g, taumax, i

    
    GLOBAL e, taumin, vhalf, s1, s2
}

PARAMETER {
    gmax= 0.001 (S/cm2)
    e = -35 (mV)

    
    
    vhalf = -78 (mV)
    s1 = -13 (mV)
    s2 = 14 (mV)
    taumax = 1350 (1)
    taumin = 10 (ms)
}

ASSIGNED { 
    v (mV)

    i (mA/cm2)
    ninf
    tau (ms)
    g (S/cm2)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g  = gmax*n
    i  = gmax*n*(v-e)
}



INITIAL {
    settables(v)
    n = ninf
}

DERIVATIVE states { 
    settables(v)      
    n' = (ninf - n)/(taumax*tau+taumin)
}


PROCEDURE settables(v (mV)) {
    
    

    TABLE ninf, tau DEPEND vhalf, s1, s2
          FROM -200 TO 50 WITH 750

    
    
    
    ninf = 1/(1+(exp((vhalf-v)/s1)))

    
    
    
    
    
    tau = (4 (ms))/(1+exp((vhalf-v)/s2))*ninf

}