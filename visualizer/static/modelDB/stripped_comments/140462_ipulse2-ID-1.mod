NEURON {
	POINT_PROCESS Ipulse2
	RANGE del, dur, per, num, amp, i
	ELECTRODE_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del (ms)
	dur (ms) <0, 1e9>	
	per (ms) <0, 1e9>	
	num			
	amp (nA)		
}

ASSIGNED {
	ival (nA)
	i (nA)
	on
	tally			
}

INITIAL {
	if (dur >= per) {
		per = dur + 1 (ms)
		printf("per must be longer than dur\n")
UNITSOFF
		printf("per has been increased to %g ms\n", per)
UNITSON
	}
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
			
			net_send(dur, 1)
		} else {
			
			ival = 0
			on = 0
			if (tally > 0) {
				
				net_send(per - dur, 1)
				tally = tally - 1
			}
		}
	}
}