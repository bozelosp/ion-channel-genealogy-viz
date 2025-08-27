NEURON {
    POINT_PROCESS tmGabaA
    RANGE tau1, tau2, e, i, q
    RANGE tau, tauR, tauF, U, u0
    RANGE base, f_gaba
    POINTER pka
    NONSPECIFIC_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
}

PARAMETER {
    tau1= 0.5 (ms)
    tau2 = 7.5 (ms)  
    e = -60 (mV)
    tau = 3 (ms)
    tauR = 500 (ms)  
    tauF = 0 (ms)    
    U = 0.1 (1) <0, 1>
    u0 = 0 (1) <0, 1>
    q = 2
    base   = 0.0      
	f_gaba = 0.0      
}

ASSIGNED {
    v (mV)
    i (nA)
    g (uS)
    factor
    x
    pka (1)
}

STATE {
    A (uS)
    B (uS)
}

INITIAL {
    LOCAL tp
    A = 0
    B = 0
    tp = (tau1*tau2)/(tau2-tau1) * log(tau2/tau1)
    factor = -exp(-tp/tau1) + exp(-tp/tau2)
    factor = 1/factor
    tau1 = tau1/q
    tau2 = tau2/q
}

BREAKPOINT {
    SOLVE state METHOD cnexp
    g = B - A
    i = modulation(f_gaba)*g*(v - e)
}

DERIVATIVE state {
    A' = -A/tau1
    B' = -B/tau2
}

NET_RECEIVE(weight (uS), y, z, u, tsyn (ms)) {
    INITIAL {
        y = 0
        z = 0
        u = u0
        tsyn = t
    }
    z = z*exp(-(t-tsyn)/tauR)
    z = z + (y*(exp(-(t-tsyn)/tau) - exp(-(t-tsyn)/tauR)) / (tau/tauR - 1) )
    y = y*exp(-(t-tsyn)/tau)
    x = 1-y-z
    if (tauF > 0) {
        u = u*exp(-(t-tsyn)/tauF)
        u = u + U*(1-u)
    } else {
        u = U
    }
    A = A + weight*factor*x*u
    B = B + weight*factor*x*u
    y = y + x*u
    tsyn = t
}

FUNCTION modulation(modFact) {
    
    
    
    modulation = 1 + modFact * (pka - base)
    
}