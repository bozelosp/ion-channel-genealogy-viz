NEURON {
    POINT_PROCESS NMDA_dsyn
    RANGE gmax
    NONSPECIFIC_CURRENT i
    RANGE tau1,tau2, factor
    RANGE g, mgB
    RANGE t_last
    GLOBAL Deadtime, Prethresh
    GLOBAL vmin, vmax
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (nS) = (nanosiemens)
}

PARAMETER {
    gmax = 1.5 (nS)
    tau1 = 75 (ms) 
    tau2 = 2 (ms)
    Prethresh = 0 (mV)  
    Deadtime = 0 (ms)  
    vmin = -120 (mV)
    vmax = 100 (mV)
}

ASSIGNED { 
    i (nA)  
    g (nS)  
    v (mV)  
    factor
    tp (ms)
    t_last (ms)  
}

STATE {
    A  
    B  
    mgB  
}

INITIAL {
    i=0 
    A=0
    B=0
    t_last = -1 - Deadtime
    tp = (tau2*tau1)/(tau1 - tau2) * log(tau1/tau2)
    factor = -exp(-tp/tau2) + exp(-tp/tau1)
    factor = 1/factor
    mgblock(v)
}    

BREAKPOINT {  
    mgblock(v)
    SOLVE state METHOD cnexp
    g=gmax*mgB*(A-B)
    i=(1e-3)*g*v
}

DERIVATIVE state {
    A'=-A/tau1
    B'=-B/tau2
}

PROCEDURE mgblock(v(mV)) {
    TABLE mgB
    FROM vmin TO vmax WITH 400
    mgB = 1/(1 + 0.3*exp(-0.1 (/mV) *(v)))
}

NET_RECEIVE(trgr) {
    if (t_last + Deadtime <= t) {
        t_last = t
        A = A + trgr*factor
        B = B + trgr*factor
    }
}