INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS hybrid
	ELECTRODE_CURRENT i
	RANGE delay, dur, vc_amp, tot_dur, rs, vc, i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	rs = 1 (megohm) <1e-9, 1e9>
    delay = 4000 (ms)
    dur = 4000 (ms)
    vc_amp = -70 (mV)
    tot_dur = 25000 (ms)

}

ASSIGNED {
	v (mV)	
	i (nA)
	vc (mV)
	tc2 (ms)
	tc3 (ms)
	on
}

INITIAL {
	
	on = 0
    i = 0
}

BREAKPOINT {
	SOLVE icur METHOD after_cvode
	vstim()
}

PROCEDURE icur() {
	if (on) {
		i = (vc - v)/rs
	}else{
		i = 0
	}
}



PROCEDURE vstim() {
    at_time(delay)
    at_time(delay+dur)
    
	if (t < delay) {
		on = 1
        vc = vc_amp
	}else if (t < delay+dur) {
        on = 0
        vc = 0
	}else {
		vc = vc_amp
		on = 1
	}
	icur()
}