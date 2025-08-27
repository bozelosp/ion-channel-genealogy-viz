DEFINE MAXPULSES 1000		

NEURON {
	POINT_PROCESS Ipulse 
	RANGE del, dur, per, num, amp, ampli,dc, i, x,pcount ,ix
	ELECTRODE_CURRENT ix
}

UNITS {
	(nA) = (nanoamp)
    PI = (pi) (1)
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
  x
  ix
  ampli
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
  x = ampli * sin(2*PI*0.001*(t-125))
  ix=i+x
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