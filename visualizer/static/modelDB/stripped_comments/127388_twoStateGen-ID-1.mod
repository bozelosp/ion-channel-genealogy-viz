NEURON	{ 
  ARTIFICIAL_CELL twoStateGen

  RANGE start
  RANGE bMean, bStd, iMean, iStd, bDurMean, bDurStd
  RANGE p
}

PARAMETER {
	start		= 0 (ms)	

	bMean		= 5 (ms) 	
	bStd		= 0
	iMean		= 100 (ms)	
	iStd		= 0
	bDurMean	= 50 (ms)   
	bDurStd	= 0


	p		= 1		
}


ASSIGNED {
	event (ms)
	on
}


PROCEDURE seed(x) {
	set_seed(x)
}


INITIAL {
	on = 0
	if (start >= 0) {
		
		
		event = start + invl(iMean,iStd) - iMean
		
		if (event < 0) {
			event = 0
		}
		net_send(event, 5)
	}

}


FUNCTION invl(mean (ms), stddev (ms)) {
	if (mean <= 0.) {
		mean = .01 (ms) 
	}

	if (stddev > 0) {
		invl = normrand(mean, stddev)
	} else {
		invl = mean
	}
}


PROCEDURE next_invl() {
	if (on == 1) {
		event = invl(bMean, bStd)
	} else {
		event = invl(iMean, iStd)
	}
	if (event < 0) {
		event = 0
	}
}


NET_RECEIVE (w) {LOCAL tmp
	if (flag == 0) { 
		if (w > 0 && on == 0) { 
			if (scop_random() <= p) {
				if (bDurMean > 0) {
					on = 1
					next_invl()
					net_send(event, 2)
					tmp = normrand(bDurMean, bDurStd)
					if (tmp < 0) {
						tmp = bDurStd
					}
					net_send(tmp, 3)	
				}
			}

		} else {
			if (w < 0 && on == 1) { 
				on = 0
			}
		}

	} else {
	if (flag == 1) {
		net_event(t)
		next_invl()
		net_send(event,1)

	} else {
	if (flag == 2) {
		net_event(t)
		next_invl()
		if (on == 1) {
			net_send(event,2)		
		}

	} else {
	if (flag == 3) {
		on = 0				

	} else {
	if (flag == 5) { 
		net_send(0,1)			
	}

	}
	}
	}
	}
}