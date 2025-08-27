NEURON {
        POINT_PROCESS RampIClamp
        RANGE del, dur, amp1, amp2
        ELECTRODE_CURRENT i
}

UNITS {
        (nA) = (nanoamp)
             }

PARAMETER {
        del=5	(ms)
        dur=200	(ms)
        amp1=1 (nA)
		amp2=1 (nA)
}

ASSIGNED {
        i (nA)
}

BREAKPOINT {
		at_time(del)
		at_time(del + dur)

		if (t < del) {
			i=0	
		}else{  
			if (t < del+dur) {
				i=amp1+(amp2-amp1)*(t-del)/dur
			}else{  
				i=0
			}
		}
	}