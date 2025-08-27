NEURON {
    POINT_PROCESS GapPC
    RANGE r, i, vgap
    ELECTRODE_CURRENT i
}

PARAMETER {
    r = 1428 (megohm)  
}

ASSIGNED {
    v     (mV)
    vgap  (mV)
    i     (nA)
}

BREAKPOINT {
    i = (vgap - v) / r
}