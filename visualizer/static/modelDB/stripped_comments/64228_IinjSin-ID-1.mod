NEURON {
	POINT_PROCESS IinjSin
	RANGE del, ton, toff, num, amp,teFreq,ssI, i
	ELECTRODE_CURRENT i
}

UNITS {
	(pA) = (picoamp)
        (nA) = (nanoamp)
}

PARAMETER {
	del  = 100 (ms)
	ton  = 500 (ms) <0, 1e9>	
	toff = 1000 (ms) <0, 1e9>	
	num  = 5			
	amp  = 10 (pA)		
        tcFreq = 5 (Hz)   
        ssI  = 40 (pA)     
}

ASSIGNED {
        Ncount     
	ival (nA)
	i (nA)
	on
	tally			
	tr (ms)   
        ssInA (nA)
}

INITIAL {
	i = 0
	ival = 0
	tally = num
        Ncount=0
        ssInA=ssI*0.001
	if (tally > 0) {
		net_send(del, 1)
		on = 0
		tally = tally - 1
	}
}

BREAKPOINT {

        tr=t-del-(ton+toff)*(Ncount-1)
        if (on ==1) { 
	i = ssInA+ival*sin(2*3.14*tcFreq*tr/1000)
        } else {
        i = ssInA+ival
        }
         
}

NET_RECEIVE (w) {
	
	if (flag == 1) {
		if (on == 0) {
			
                        Ncount=Ncount+1
			ival = amp*0.001
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