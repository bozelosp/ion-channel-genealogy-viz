INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE NSTEP 3

NEURON {
	POINT_PROCESS SEClamp1
	ELECTRODE_CURRENT i
	RANGE dur1, amp1, dur2, amp2, dur3, amp3, dur4, amp4, rs, vc, i
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
	dur4 (ms) <0,1e9> amp4 (mV)
}

ASSIGNED {
	v (mV)	
	i (nA)
	vc (mV)
	tc2 (ms)
	tc3 (ms)
	tc4 (ms)
	on
}

INITIAL {
	tc2 = dur1 + dur2
	tc3 = tc2 + dur3
	tc4 = tc3 + dur4
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
	if (dur1) {at_time(dur1)}
	if (dur2) {at_time(tc2)}
	if (dur3) {at_time(tc3)}
	if (dur4) {at_time(tc4)}
	if (t < dur1) {
		vc = amp1
	}else if (t < tc2) {
		vc = amp2
	}else if (t < tc3) {
		vc = amp3
	}else if (t < tc4) {
		vc = amp4
	}else {
		vc = 0
		on = 0
	}
	icur()
}