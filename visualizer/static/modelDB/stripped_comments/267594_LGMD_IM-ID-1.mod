UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}


NEURON {
    THREADSAFE
    
    SUFFIX M

    USEION k READ ek WRITE ik

    
    RANGE vhalf, gmax, tau, g

    
    GLOBAL taumax, taumin, s1, s2
}

PARAMETER {
    gmax= 0.0003 (S/cm2)

    vhalf = -50 (mV)	
    s1 = 15		(mV)	
	s2 = -20	(mV)
	taumax = 80 (ms)
	taumin = 15 (ms)
	aop = 0		(1)	< 0, 1 >	
}

ASSIGNED { 
    v	(mV)
    ek	(mV)
    
    ik	(mA/cm2)
    ninf	(1)
    tau (ms)
    g 	(S/cm2)
}

STATE {
    n
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g  = gmax*(n*(1-aop)+aop)
    ik  = g*(v-ek)
}



INITIAL {
    settables(vhalf-v)
    n = ninf
}

DERIVATIVE states { 
    settables(vhalf-v)      
    n' = (ninf - n)/tau
}


PROCEDURE settables(df (mV)) {
	
    TABLE ninf, tau DEPEND s1, s2, taumax, taumin
          FROM -100 TO 100 WITH 600

    
    ninf = 1/(1+(exp((df)/s1)))

    
    
    tau = 4*(taumax-taumin)/(1+exp((df)/s2))*ninf+taumin
}