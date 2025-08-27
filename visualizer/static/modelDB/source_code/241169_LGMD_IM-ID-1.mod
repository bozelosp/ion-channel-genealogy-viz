TITLE IM channel for LGMD


UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

NEURON {
    THREADSAFE
    : note - every variable accessible in NEURON will have the suffix _M
    SUFFIX M

    USEION k READ ek WRITE ik

    : these variables can be accessed as compartment.rangevar_M
    RANGE vhalf, gmax, tau, g

    : these will be accessed as taumax_M, vhalf_M, s1_M, s2_M  
    GLOBAL taumax, taumin, s1, s2
}

PARAMETER {
    gmax= 0.0003 (S/cm2)

    vhalf = -48	(mV)
    s1 = 12		(mV)
    s2 = -8		(mV)
    taumax = 20	(ms)
    taumin =2.5	(ms)
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
    g  = gmax*n
    ik  = g*(v-ek)
}

: calls the function settables below, then 
: set the steady state value of IM activation 
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
          FROM -60 TO 60 WITH 600

    : steady-state activation of IM in mV
    ninf = 1/(1+(exp((df)/s1)))

    : steady-state IM time constant
    tau = 4*(taumax-taumin)/(1+exp((df)/s2))*ninf+taumin
}


