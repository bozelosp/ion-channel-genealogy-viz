NEURON {
	POINT_PROCESS TriggeredIClamp
	RANGE dur, amp, i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	dur (ms)	<0,1e9>
	amp (nA)
}

ASSIGNED { i (nA) ilocal (nA)}

INITIAL {
	i = 0
	ilocal = 0
}

BREAKPOINT {
		i = ilocal
}

NET_RECEIVE(w) {
	if (flag == 0) {
		ilocal = ilocal + amp
		net_send(dur, 1)
	}else{
		ilocal = ilocal - amp
	}
}