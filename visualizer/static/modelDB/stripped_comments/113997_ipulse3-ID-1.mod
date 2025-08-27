DEFINE MAXPULSES 1000		

NEURON {
	POINT_PROCESS Ipulse3
	RANGE del, dur, per, num, amp, dc, i, pcount
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
	dc (nA)			
}

ASSIGNED {
	amp[MAXPULSES]		
	ival (nA)
	i (nA)
	on
	tally			
	pcount			
}

INITIAL {
	pcount = 0
	if (dur >= per) {
		per = dur + 1 (ms)
		printf("per must be longer than dur\n")
UNITSOFF
		printf("per has been increased to %g ms\n", per)
UNITSON
	}
	i = 0
	ival = dc
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
			
			ival = amp[pcount]+dc
			pcount = pcount + 1
			on = 1
			
			net_send(dur, 1)
		} else {
			
			ival = dc
			on = 0
			if (tally > 0) {
				
				net_send(per - dur, 1)
				tally = tally - 1
			}
		}
	}
}