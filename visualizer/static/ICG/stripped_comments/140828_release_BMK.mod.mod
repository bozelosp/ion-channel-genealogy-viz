INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX rel
	RANGE T, del, dur, amp
}

UNITS {
	(mM) = (milli/liter)
}

PARAMETER {
	del (ms)
	dur (ms)	<0,1e9>
	amp (mM)
}

ASSIGNED { T (mM)
}


INITIAL {
	T = 0
}

BREAKPOINT {
	at_time(del)
	at_time(del+dur)

	if (t < del + dur && t > del) {
		T = amp
	}else{
		T = 0
	}
}