UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}



NEURON {
    THREADSAFE
    
    SUFFIX M

    USEION k READ ek WRITE ik

    
    
    
    RANGE gmax, g

    
    GLOBAL taumax, taumin, vhalf, s1, s2
}

PARAMETER {

    gmax= 0.001 (S/cm2)
    vhalf = -51	(mV)
    v2 = -60	(mV)
    s1 = 10.5	(mV)
    s2 = -11 	(mV)
    taumax = 60 (ms)
    taumin = 2	(ms)
}

ASSIGNED { 
    v (mV)
    ek (mV)
    
    ik (mA/cm2)
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
    ik = g*(v-ek)
}



INITIAL {
    settables(v)
    n = ninf
}

DERIVATIVE states { 
    settables(v)      
    n' = (ninf - n)/tau
}


PROCEDURE settables(v (mV)) {
    
    

    TABLE ninf, tau DEPEND vhalf, s1, s2, taumax, taumin
          FROM -100 TO 50 WITH 750

    
    ninf = 1/(1+(exp((vhalf-v)/s1)))

    
    
	
	tau = 2*taumax/( exp((vhalf-v)/s2) + exp((vhalf-v)/s1) ) + taumin
}