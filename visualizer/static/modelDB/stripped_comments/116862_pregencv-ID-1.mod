INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON	{ 
  POINT_PROCESS SpikeGenerator
  RANGE x, spk, lastspk
  RANGE fast_invl, slow_invl, burst_len, start, end
  RANGE noise
  GLOBAL dummy 
}

PARAMETER {
	fast_invl	= 1		
	slow_invl	= 50		


	burst_len	= 10		
	start		= 50		
	end		= 1e10		
	noise		= 0		
}

ASSIGNED {
	x
	burst
	event
	burst_off
	burst_on
	toff
	dummy
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	toff = 1e9
	x = -90
	burst = 0
	event = start
	event_time()
	generate()
}	

BREAKPOINT {
	SOLVE generate METHOD cvode_t
}

FUNCTION interval(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 
	}
	if (noise == 0) {
		interval = mean
	}else{
		interval = (1. - noise)*mean + noise*exprand(mean)
	}
}

PROCEDURE event_time() {
	if (burst != 0.) {
		event = event + interval(fast_invl)
		if (event > burst_on + burst_off) {
			burst = 0.
		}
	}else{
		burst = 1.



		event = event + interval(slow_invl)
		burst_on = event
		burst_off = interval((burst_len - 1)*fast_invl)-1e-6
	}
	if (event > end) {
		event = -1e5
	}
}

PROCEDURE generate() {
	if (at_time(event)) {
		x = 20
		toff = event + .1
		event_time()
	}
	if (at_time(toff)) {
		x = -90
	}
}