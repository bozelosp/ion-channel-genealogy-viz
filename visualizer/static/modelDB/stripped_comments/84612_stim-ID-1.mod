INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS MyIClamp
	RANGE del,dur,amp,amp2,i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del  (ms)
	dur  (ms)	<0,1e9>
	amp  (nA)
	amp2 (nA)
}
ASSIGNED { i (nA) }

INITIAL {
	i = amp
}

BREAKPOINT {
	at_time(del)
	at_time(del+dur)

	if (t < del + dur && t > del) {
		i = amp +amp2
	}else{
		i = amp
	}
}