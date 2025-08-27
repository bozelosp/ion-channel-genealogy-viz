INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 3

NEURON {
        POINT_PROCESS SEVClamp
        ELECTRODE_CURRENT i
        RANGE dur, amp, rs, vc, i
}

UNITS {
        (nA) = (nanoamp)
        (mV) = (millivolt)
        (uS) = (micromho)
}


PARAMETER {
        v (mV)
        rs = 1 (megohm)		
}

ASSIGNED {
        i (nA)
        vc (mV)
        ic (nA)
        tc2 (ms)
        tc3 (ms)
	dur[NSTEP] (ms)
	amp[NSTEP] (mV)
        on
}

INITIAL {
        tc2 = dur[0] + dur[1]
        tc3 = tc2 + dur[2]
        on = 0
}

BREAKPOINT {
        SOLVE vstim
        if (on) {
                i = (vc - v)/rs
        }else{
                i = 0
        }
}

PROCEDURE vstim() {
        on = 1
        if (t < dur[0]) {
                vc = amp[0]
        }else if (t < tc2) {
                vc = amp[1]
        }else if (t < tc3) {
                vc = amp[2]
        }else {
                vc = 0
                on = 0
        }
        if (on) {
        }else{
                ic = 0
        }
        VERBATIM
        return 0;
        ENDVERBATIM
}