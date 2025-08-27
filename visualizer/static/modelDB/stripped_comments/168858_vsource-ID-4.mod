NEURON {
	POINT_PROCESS Vsource
	ELECTRODE_CURRENT i
	RANGE toff, amp, rs, vc, i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}


PARAMETER {
	rs = 1 (megohm) <1e-9, 1e9>
	toff (ms) 	  amp (mV)
}

ASSIGNED {
	v (mV)	
	i (nA)
	vc (mV)
	on
}

INITIAL {
	on = 1
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
	if (toff) {at_time(toff)}
	if (t < toff) {
		vc = amp
	}else {
		vc = 0
		on = 0
	}
	icur()
}