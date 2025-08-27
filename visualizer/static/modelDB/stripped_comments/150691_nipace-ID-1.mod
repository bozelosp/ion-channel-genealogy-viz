NEURON {
	POINT_PROCESS NIClamp
	RANGE del, dur, amp, del1, n, i
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del = 100 (ms)	<0,1e9>
	dur = .1 (ms)	<0,1e9>
	amp = 2 (nA)
	del1 = 1000 (ms) <1e-9,1e9> 
	n = 100 	<0,1e9> 
}
ASSIGNED { i (nA) cnt a(nA)}

INITIAL {
	i = 0
	cnt = 0
	net_send(del, 1)
}

BREAKPOINT {
	i = a
}

NET_RECEIVE(w) {
	if (flag == 1) {
		a = amp
		cnt = cnt + 1
		net_send(dur, 2)
	}else if (flag == 2) {
		a = 0
		if (cnt < n) {
			net_send(del1, 1)
		}
	}
}