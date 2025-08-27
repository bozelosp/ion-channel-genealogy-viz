NEURON {
	POINT_PROCESS IinjLT
	RANGE del, ton, toff, num, amp,ssI,i
	ELECTRODE_CURRENT i
}

UNITS {
	(pA) = (picoamp)
    (nA) = (nanoamp)
}

PARAMETER {
	del  = 500 (ms)
	ton  = 8000 (ms) <0, 1e9>	
	toff = 0 (ms) <0, 1e9>	
	num  = 1			
	amp  = 0 (pA)	  
    ssI  = 40 (pA)     
        
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
			Part1=32*( 1-exp(- (tr/1000 )/(0.05/1.2*ton)  ) )
			Part2=-33/(   1+exp(-   (  (tr/(1000)) -((3.8/1.2)*ton)   )/((0.45/1.2)*ton)    ) )
			Part3=1-exp(  - (tr/(1000))  /((0.8/1.2)*ton) )
			i=ssInA-ival*(Part1+Part2+Part3)/33.0
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