NEURON {
	POINT_PROCESS Pulses
	RANGE del, dur, amp, i, npulses, period, n
	ELECTRODE_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del (ms)
	dur (ms)	<0,1e9>
	amp (nA)
	npulses
	period (ms)
}

ASSIGNED { i (nA)
		   n
		 }

INITIAL {
	i = 0
	n=0
}

BREAKPOINT {
	
	if (n==0) {
		if((t >= del) && (t<(del+dur))) { 
			i = amp
		} 
		if(t >= (del+dur)){ 
			if(i==amp){
				n=n+1
			}
			i=0
		}	
	} else if (n<npulses) {
		if(t>=(del+n*period) && t<(del+n*period+dur)) {
			i=amp
		}
		if(t >= (del+n*period+dur)){ 
			if(i==amp){
				n=n+1
			}
			i=0
		}	
	}
}