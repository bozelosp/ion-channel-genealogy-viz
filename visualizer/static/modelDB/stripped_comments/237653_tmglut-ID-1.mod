NEURON {
    THREADSAFE
    POINT_PROCESS tmGlut
    RANGE tau1_ampa, tau2_ampa, tau1_nmda, tau2_nmda
    RANGE g_ampa, g_nmda, i_ampa, i_nmda, nmda_ratio
    RANGE e, g, i, q, mg
    RANGE tau, tauR, tauF, U, u0
    RANGE ca_ratio_ampa, ca_ratio_nmda
    RANGE base, f_ampa, f_nmda
    POINTER pka
    NONSPECIFIC_CURRENT i
    USEION cal WRITE ical VALENCE 2
}

UNITS {
    (nA) = (nanoamp)
    (mV) = (millivolt)
    (uS) = (microsiemens)
    (mM) = (milli/liter)
}

PARAMETER {
    tau1_ampa= 2.2 (ms)
    tau2_ampa = 11.5 (ms)  
    tau1_nmda= 5.63 (ms)
    tau2_nmda = 320 (ms)  
    nmda_ratio = 0.5 (1)
    e = 0 (mV)
    tau = 3 (ms)
    tauR = 100 (ms)  
    tauF = 0 (ms)  
    U = 0.3 (1) <0, 1>
    u0 = 0 (1) <0, 1>
    ca_ratio_ampa = 0.005
    ca_ratio_nmda = 0.1
    mg = 1 (mM)
    q = 2
    base   = 0.0      
	f_ampa = 0.0      
	f_nmda = 0.0      
}

ASSIGNED {
    v (mV)
    i (nA)
    i_ampa (nA)
    i_nmda (nA)
    ical (nA)
    ical_ampa (nA)
    ical_nmda (nA)
    g (uS)
    g_ampa (uS)
    g_nmda (uS)
    factor_ampa
    factor_nmda
    x
    pka (1)
}

STATE {
    A_ampa (uS)
    B_ampa (uS)
    A_nmda (uS)
    B_nmda (uS)
}

INITIAL {
    LOCAL tp_ampa, tp_nmda
    A_ampa = 0
    B_ampa = 0
    tp_ampa = (tau1_ampa*tau2_ampa)/(tau2_ampa-tau1_ampa) * log(tau2_ampa/tau1_ampa)
    factor_ampa = -exp(-tp_ampa/tau1_ampa) + exp(-tp_ampa/tau2_ampa)
    factor_ampa = 1/factor_ampa
    tau1_ampa = tau1_ampa/q
    tau2_ampa = tau2_ampa/q
    A_nmda = 0
    B_nmda = 0
    tp_nmda = (tau1_nmda*tau2_nmda)/(tau2_nmda-tau1_nmda) * log(tau2_nmda/tau1_nmda)
    factor_nmda = -exp(-tp_nmda/tau1_nmda) + exp(-tp_nmda/tau2_nmda)
    factor_nmda = 1/factor_nmda
    tau1_nmda = tau1_nmda/q
    tau2_nmda = tau2_nmda/q
}

BREAKPOINT {
    LOCAL itotal, mggate
    SOLVE state METHOD cnexp
    mggate = 1 / (1 + exp(-0.062 (/mV) * v) * (mg / 3.57 (mM)))
    g_ampa = B_ampa - A_ampa
    itotal = modulation(f_ampa) * g_ampa*(v - e)
    ical_ampa = ca_ratio_ampa*itotal
    i_ampa = itotal - ical_ampa
    g_nmda = B_nmda - A_nmda
    itotal = modulation(f_nmda) * g_nmda*(v - e)*mggate
    i_nmda = itotal - ical_nmda
    ical = ical_ampa + ical_nmda
    i = i_ampa + i_nmda
    g = g_ampa + g_nmda
}

DERIVATIVE state {
    A_ampa' = -A_ampa/tau1_ampa
    B_ampa' = -B_ampa/tau2_ampa
    A_nmda' = -A_nmda/tau1_nmda
    B_nmda' = -B_nmda/tau2_nmda
}

NET_RECEIVE(weight (uS), y, z, u, tsyn (ms)) {
    LOCAL weight_ampa, weight_nmda
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
    weight_ampa = weight
    weight_nmda = weight_ampa*nmda_ratio
    A_ampa = A_ampa + weight_ampa*factor_ampa*x*u
    B_ampa = B_ampa + weight_ampa*factor_ampa*x*u
    A_nmda = A_nmda + weight_nmda*factor_nmda*x*u
    B_nmda = B_nmda + weight_nmda*factor_nmda*x*u
    y = y + x*u
    tsyn = t
}

FUNCTION modulation(factor) {
    
    
    
    modulation = 1 + factor * (pka - base)
    
}