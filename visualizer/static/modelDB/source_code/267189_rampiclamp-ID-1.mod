COMMENT
Since this is an electrode current, positive values of i depolarize the cell
and in the presence of the extracellular mechanism there will be a change
in vext since i is not a transmembrane current but a current injected
directly to the inside of the cell.
By Kelvin and Ted, 2007
https://www.neuron.yale.edu/phpBB/viewtopic.php?t=897
ENDCOMMENT

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