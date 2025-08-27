NEURON {
	POINT_PROCESS Gpulse1
	RANGE del, ton, toff, num, e, gmax, g, i
	NONSPECIFIC_CURRENT i
}

UNITS {
	(mV) = (millivolt)
	(nA) = (nanoamp)
	(uS) = (microsiemens)
}

PARAMETER {
	del (ms)
	ton (ms) <0, 1e9>	
	toff (ms) <0, 1e9>	
	num			
	e (mV)			
	gmax (uS) <0, 1e9>	
}

ASSIGNED {
	v (mV)
	g (uS)		

	i (nA)
	on
	tally		
}

INITIAL {
	g = 0
	i = 0

	tally = num
	if (tally > 0) {
		net_send(del, 1)
		on = 0
		tally = tally - 1
	}
}

BREAKPOINT {

	i = g*(v-e)
}

NET_RECEIVE (w) {
	
	if (flag == 1) {
		if (on == 0) {
			

			g = gmax
			on = 1
			
			net_send(ton, 1)
		} else {
			

			g = 0
			on = 0
			if (tally > 0) {
				
				net_send(toff, 1)
				tally = tally - 1
			}
		}
	}
}