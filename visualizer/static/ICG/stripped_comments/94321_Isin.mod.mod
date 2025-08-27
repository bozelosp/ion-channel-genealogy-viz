NEURON {
    POINT_PROCESS Isin
    RANGE del, dur, amp, i, freq, off
    ELECTRODE_CURRENT i
}
UNITS {
    (nA) = (nanoamp)
}

PARAMETER {
    off (nA)
    del (ms)
    dur (ms)    <0,1e9>
    amp (nA)
    freq (Hz)
}
ASSIGNED { i (nA) }

INITIAL {
    i = 0
}

BREAKPOINT {
    at_time(del)
    at_time(del+dur)

    if (t < del + dur && t >= del) {
        i = off + amp * sin(0.0062831853071795866*freq*t)
    }else{
        i = off
    }
}