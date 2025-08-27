NEURON {
    POINT_PROCESS Synapse_HorPR
    RANGE V_pre, v_th, v_slope, g_max, i
    NONSPECIFIC_CURRENT i
}

PARAMETER {
    v_th    = -30.88  (mV)
    v_slope = 5	(mV)
    g_max   = 0.250   (umho)
    
}

ASSIGNED {
    V_pre  (mV)
    i      (nA)
    g      (nA)
}

BREAKPOINT {
    g = tanh( (V_pre - v_th) / v_slope ) + 1
    i = -g_max * 0.5 * g
}