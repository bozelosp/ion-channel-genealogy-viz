TITLE Ih channel for LGMD
: this channel is identical to that implemented in Matlab code
: based on Richard's current and voltage clamp data (Sept 13)
: this is a second version were the state variables are ninf and tau 
: instead of nalpha and nbeta to allow for inspection of these variables 
: during simulations
: Modified on 04/08/14 to correct the apparent mix-up between tau and tau_max
: tau_max was RANGE and tau GLOBAL when it should be the opposite. 

UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
    (S) = (siemens)
}

: gmax and g are range variables (i.e., can change in different compartments
: while e is global
NEURON {
    THREADSAFE
    SUFFIX h
    NONSPECIFIC_CURRENT i
    RANGE gmax, tau, g, taumax, i
    GLOBAL e, taumin, vhalf, s1, s2
}

PARAMETER {
    gmax= 0.001 (S/cm2)
    e = -37		(mV)
    vhalf = -77 (mV)
    s1 = -12	(mV)
    s2 = 13		(mV)
    taumax =1350 (ms)
    taumin = 10 (ms)
}

ASSIGNED { 
    v	(mV)
    i	(mA/cm2)
    ninf	(1)
    tau	(1)
    g	(S/cm2)
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
    :LOCAL vhalf, s1, s2, taumax
    :local variables take units of right hand side, see below

    TABLE ninf, tau DEPEND vhalf, s1, s2
          FROM -120 TO 20 WITH 500

    : steady-state activation of Ih in mV
    ninf = 1/(1+(exp((vhalf-v)/s1)))

    : steady-state Ih time constant
    : slope in mV and time constant in ms
    tau = 4/(1+exp((vhalf-v)/s2))*ninf

}


