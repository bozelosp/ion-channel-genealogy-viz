NEURON {
	POINT_PROCESS IinjLT_cone
	RANGE del, ton, toff, num, amp,ssI,i
	ELECTRODE_CURRENT i
}

UNITS {
	(pA) = (picoamp)
        (nA) = (nanoamp)
}

PARAMETER {
	del  = 1000 (ms)
	ton  = 8000 (ms) <0, 1e9>	
	toff = 1000 (ms) <0, 1e9>	
	num  = 1			
	amp  = 20 (pA)	  
        ssI  = 20 (pA)     
        
}

ASSIGNED {
        Ncount     
	ival (nA)
	i (nA)
	on
	tally			
	tr (ms)   
        Part1
        Part2
        Part3
	Part4
        ssInA (nA)
        ampnA (nA)
}

INITIAL {
	i = 0
        ssInA=0.001*ssI
        ampnA=amp*0.001
	ival = 0
	tally = num
        Ncount=0
	if (tally > 0) {
		net_send(del, 1)
		on = 0
		tally = tally - 1
	}
}

BREAKPOINT {

        tr=t-del-(ton+toff)*(Ncount-1)
        if (on ==1) { 
	 
          Part1=(1.01)*( 1-exp(-tr/60))^3
          Part2=-(0.07)/(1+exp(-(tr-555)/120) )
          Part3=-0.7/(1+exp(-(tr-1500)/200))^1.3
	  Part4= -(0.24)/((1+exp(-(tr-1200)/80)))^2
          i=ssInA-ival*(Part1+Part2+Part3+Part4)
                   
        } else {
        i = ssInA+ival
        }
         
}

NET_RECEIVE (w) {
	
	if (flag == 1) {
		if (on == 0) {
			
                        Ncount=Ncount+1
			ival = ampnA
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