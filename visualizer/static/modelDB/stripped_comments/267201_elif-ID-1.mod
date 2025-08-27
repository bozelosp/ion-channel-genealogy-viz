NEURON {
    POINT_PROCESS eLIF
    RANGE V_reset, t_ref, V_th, spikewidth, V_peak
    RANGE epsilon
    RANGE g_L, E_0, E_u, E_d, E_f
    RANGE tau_e, epsilon_0, epsilon_c, delta, alpha
    NONSPECIFIC_CURRENT i, I_e
}

UNITS {
    (mV) = (millivolt)
    (pA) = (picoamp)
    (nS) = (nanosiemens)
}

PARAMETER {
    V_th = -50   (mV)   
    V_reset  = -60   (mV)   
    V_peak  = 0.   (mV)   
    t_ref = 2     (ms)   
    gon     = 1e9   (nS)   
    spikewidth = 1e-12 (ms) 

    E_0 = -55. (mV)
    E_u = -50. (mV)
    E_d = -35. (mV)
    E_f = -45. (mV)

    epsilon_0 = 0.5
    epsilon_c = 0.15
    alpha = 1. <0, 1e9>
    delta = 0.02
    tau_e = 500. (ms) <0, 1e9>

    g_L	    = 30.   (nS)

    I_e = 0 (pA)
}


ASSIGNED {
    i (pA)
    irefrac (pA)
    grefrac (nS)
    refractory
}

STATE {
    w  (pA)
    epsilon
}

INITIAL {
    grefrac = 0
    net_send(0,4)
    epsilon = alpha*epsilon_0
}

FUNCTION E_L (epsilon) (mV) {
    E_L = E_0 + (E_u - E_0)*(1 - epsilon/epsilon_0)
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit  
    irefrac = grefrac*(v-V_reset)
    i = ( g_L*(v - E_L(epsilon)) + irefrac - I_e)
}

DERIVATIVE states {		
    epsilon' = ((1-epsilon/(alpha*epsilon_0))*(1-epsilon/(alpha*epsilon_0))*(1-epsilon/(alpha*epsilon_0)) - (v-E_f)/(E_d-E_f)) / tau_e
}

FUNCTION threshcrossing (v (mV), epsilon) {
    if ((v > V_th) && (epsilon > epsilon_c)) {
        threshcrossing = 1
    }
    else {
        threshcrossing = -1
    }
}

NET_RECEIVE (weight) {
    if (flag == 1) {        
        v = V_peak
        epsilon = epsilon - delta
        net_send(spikewidth, 2)
        net_event(t)
        
    } else if (flag == 2) { 
        v = V_reset
        grefrac = gon
        if (t_ref > spikewidth) {
            net_send(t_ref-spikewidth, 3)
        } else { 
            grefrac = 0
        }
        
    } else if (flag == 3) { 
        v = V_reset
        grefrac = 0
        
    } else if (flag == 4) { 
        WATCH ( threshcrossing(v, epsilon) > 0 ) 1
    }
}