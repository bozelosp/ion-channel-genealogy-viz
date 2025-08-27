NEURON {
    POINT_PROCESS mAdExp
    RANGE V_reset, t_ref, V_spike, V_th, V_peak, spikewidth, C_m
    RANGE w, winit, epsilon
    RANGE a, b, tau_w, g_L, Delta_T, E_0, E_u, E_d, E_f
    RANGE gamma, tau_e, epsilon_0, epsilon_c, delta, I_KATP, alpha
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
    V_spike  = -43   (mV)   
    V_peak   = 0     (mV)   
    t_ref = 1     (ms)   
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
    gamma = 200. (pA) <0, 1e9>
    I_KATP = 1 (pA)

    a 	    = 0.  (nS)   
    b	    = 80.5 (pA)   
    tau_w   = 144    (ms)   

    g_L	    = 30.   (nS)   
    Delta_T   = 2      (mV)   

    I_e = 0 (pA)
    winit  = 0      (pA)
}


ASSIGNED {
    i (pA)
    irefrac (pA)
    iexp (pA)
    grefrac (nS)
    refractory
    spike_threshold (mV)
}

STATE {
    w  (pA)
    epsilon
}

INITIAL {
    grefrac = 0
    net_send(0,4)
    w = winit
    epsilon = alpha*epsilon_0
    if (Delta_T == 0) {
        spike_threshold = V_th
    } else {
        spike_threshold = V_spike
    }
}

FUNCTION E_L (epsilon) (mV) {
    E_L = E_0 + (E_u - E_0)*(1 - epsilon/epsilon_0)
}

BREAKPOINT {
    SOLVE states METHOD derivimplicit  
    irefrac = grefrac*(v-V_reset)
    iexp = exp_current(v, epsilon)
    i = ( g_L*(v - E_L(epsilon)) + iexp + w + irefrac - I_e)
}

DERIVATIVE states {		
    w' = (a*(v-E_L(epsilon)) - w + I_KATP*epsilon_c/(2*epsilon + epsilon_c)) / tau_w
    epsilon' = ((1-epsilon/(alpha*epsilon_0))*(1-epsilon/(alpha*epsilon_0))*(1-epsilon/(alpha*epsilon_0)) - (v-E_f)/(E_d-E_f) - w/gamma) / tau_e
}

FUNCTION exp_current(v, epsilon) {  
    if (Delta_T == 0) {
        exp_current = 0
    } else if ((v - V_th)/Delta_T > 100) {
        exp_current = -g_L*Delta_T*(epsilon-epsilon_c)*exp(99)/epsilon_0
    } else {
        exp_current = -g_L*Delta_T*(epsilon-epsilon_c)*exp((v-V_th)/Delta_T)/epsilon_0
    }
}

FUNCTION threshcrossing (v (mV), epsilon) {
    if ((v > V_spike) && (epsilon > epsilon_c)) {
        
        threshcrossing = 1
    }
    else {
        threshcrossing = -1
    }
}

NET_RECEIVE (weight) {
    if (flag == 1) {        
        v = V_peak
        w = w + b
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