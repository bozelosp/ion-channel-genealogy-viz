INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 3

NEURON {
	POINT_PROCESS jSEClamp
	ELECTRODE_CURRENT i
	RANGE dur1, amp1, dur2, amp2, dur3, amp3, rs, vc, i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	rs = 1 (megohm) <1e-9, 1e9>
	dur1 (ms) 	  amp1 (mV)
	dur2 (ms) <0,1e9> amp2 (mV)
	dur3 (ms) <0,1e9> amp3 (mV)
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
	tc2 = dur1 + dur2
	tc3 = tc2 + dur3
	on = 0
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
	on = 1
        if (t < 100) {
	    vc = amp1
	} else if (t < 5100){
	    vc = amp2 + ((t-100)*(amp3-amp2)/5000.0)
	} else if (t < 7100){
	    vc = amp1
	} else if (t < 12100){
	    vc = amp2 + ((t-7100)*(amp3-amp2)/5000.0)
	}else {
	    vc = amp1
	}
	icur()
}