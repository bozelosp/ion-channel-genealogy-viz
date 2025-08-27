NEURON {
	POINT_PROCESS Istim
	RANGE del, ton, toff, num, amp, i
	ELECTRODE_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del=1010 (ms)
	ton=1 (ms) <0, 1e9>	
	toff (ms) <0, 1e9>	
	num			
	amp=1 (nA)		
}

ASSIGNED {
	ival (nA)
	i (nA)
	on
	tally			
}

INITIAL {
	i = 0
	ival = 0
	tally = num
	if (tally > 0) {
		net_send(del, 1)
		on = 0
		tally = tally - 1
	}
}

BREAKPOINT {

	i = ival
}

NET_RECEIVE (w) {
	
	if (flag == 1) {
		if (on == 0) {
			
			ival = amp
			on = 1
			
			net_send(ton, 1)
		} else {
			
			ival = 0
			on = 0
			if (tally > 0) {
				
				net_send(toff, 1)
				tally = tally - 1
			}
		}
	}
}