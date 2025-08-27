NEURON	{ 
  POINT_PROCESS SpikeGenerator
  RANGE x, spk, lastspk
  RANGE fast_invl, slow_invl, burst_len, start, end
  RANGE noise
  GLOBAL dummy 
}

PARAMETER {
	fast_invl	= 1 (ms)	
	slow_invl	= 200 (ms)	


	burst_len	= 1		
	start		= 100 (ms)	
	end		= 1e10 (ms)	
	noise		= 0		
}

ASSIGNED {
	x
	burst
	event (ms)
	burst_off (ms)
	burst_on (ms)
	toff (ms)
	dummy
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	toff = 1e9
	x = -90
	burst = 0
	event = start - slow_invl
	event_time()
	while (event < 0) {
		event_time()
	}
	generate()
}	

BREAKPOINT {
	SOLVE generate METHOD cvode_t
}

FUNCTION interval(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) 
		
		
	}
	if (noise == 0) {
		interval = mean
	}else{
		interval = (1. - noise)*mean + noise*mean*exprand(1)
	}
}

PROCEDURE event_time() {
	if (slow_invl == 0 || (burst != 0. && burst_len > 1)) {
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