NEURON {
	POINT_PROCESS dbsStim
	RANGE del, dur, amp, pw, period, actPrb, adLat, active, i
	ELECTRODE_CURRENT i
}


UNITS {
	(nA) = (nanoamp)
}


PARAMETER {
	del = 0 (ms)
	dur = 1e9 (ms)	<0,1e9>
	amp = 0 (nA)
	pw = .1  (ms)
	period = 500 (ms)
	actPrb = 1
	adLat = 2
	active = 0
}


ASSIGNED { 
	i (nA) 
	offtime (ms)
}


INITIAL {
	offtime = period - pw
	i = 0
	if (del <= 0) {
		del = 0
	}
	if (active) {
		net_send(del, 1)		
		net_send(del+dur, 4)	
	}
}


BREAKPOINT {
	i = i
}


NET_RECEIVE (w) {LOCAL tmp
	if (flag == 0 && w >= 0) {
		
		active = 1
		net_send(del, 1)
		net_send(del+dur, 4)

	} else {
	if (flag == 1 && active) {		

		tmp = scop_random()
		if (tmp <= actPrb) {

			
			
			

			net_send(adLat, 2)		
			net_send(adLat+pw, 3)		
			net_event(t)			
		}
		net_send(adLat+offtime, 1)	

	} else {
	if (flag == 2) {
		i = amp			
	} else {
	if (flag == 3) {
		i = 0				
	} else {
	if (flag == 4) {
		active = 0
	}

	}
	}
	}
	}


}