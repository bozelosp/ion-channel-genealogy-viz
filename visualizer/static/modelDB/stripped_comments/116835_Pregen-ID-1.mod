NEURON	{ 
  POINT_PROCESS SpikeGenerator
  RANGE y
  RANGE fast_invl, slow_invl, burst_len, start, end,delay
  RANGE noise
}

PARAMETER {
	fast_invl	= 10 (ms)	
	slow_invl	= 0 (ms)	


	burst_len	= 10		
	start		= 50 (ms)	
	end		= 1e10 (ms)	
	noise		= 0		
	delay		= 4
}

ASSIGNED {
	y
	burst
	event (ms)
	burst_off (ms)
	burst_on (ms)
	toff (ms)
	on
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 1
	toff = 1e9
	y = -90
	burst = 0
	event = start - slow_invl
	
	event_time()
	while (on == 1 && event < 0) {
		event_time()
	}
	if (on == 1) {
		net_send(event, 1)
	}
}	

FUNCTION interval(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) 
		
		
	}
	if (noise == 0) {
		interval = mean
	}else{
		interval = (1. - noise)*mean + noise*(mean*exprand(1)+delay) 
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
		on = 0
	}
}

NET_RECEIVE (w) {

	if (flag == 1 && on == 1) {
		y = 20
		net_event(t)
		event_time()
		net_send(event - t, 1)
		net_send(.1, 2)
	}
	if (flag == 2) {
		y = -90
	}
}