NEURON {
    POINT_PROCESS PRInput
    RANGE delSteady, durSteady, ampSteady, i
    ELECTRODE_CURRENT i
}

UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    delSteady (ms)
    durSteady (ms)  < 0, 1e9 >
    ampSteady (nA)
}

ASSIGNED {
    i (nA)
}

BREAKPOINT {
    
    at_time(delSteady)
    at_time(delSteady + durSteady)

    if (t > delSteady && t < delSteady + durSteady) {
        i = ampSteady
    } else {
        i = 0
    }
}