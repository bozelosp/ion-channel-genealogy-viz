NEURON	{ 
  ARTIFICIAL_CELL NetStim_bgk
  RANGE interval, number, start
  RANGE noise
}

PARAMETER {
	interval	= 10 (ms) <1e-9,1e9>
	number	= 10 <0,1e9>	
	start		= 50 (ms)	
	noise		= 0 <0,1>	
}

ASSIGNED {
	event (ms)
	on
	ispike
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0
	ispike = 0
	if (noise < 0) {
		noise = 0
	}
	if (noise > 1) {
		noise = 1
	}
	if (start >= 0 && number > 0) {
		
		
		event = start + invl(interval) - interval*(1. - noise)
		
		if (event < 0) {
			event = 0
		}
		net_send(event, 3)
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = 0
		ispike = 0
	}
}

FUNCTION invl(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) 
	}
	if (noise == 0) {
		invl = mean
	}else{
		invl = (1. - noise)*mean + noise*mean*exprand(1)
	}
}

PROCEDURE next_invl() {
	if (number > 0) {
		event = invl(interval)
	}
	if (ispike >= number) {
		on = 0
	}
}

NET_RECEIVE (w) {
	if (flag == 0) { 
		if (w > 0 && on == 0) { 
			init_sequence(t)
			
			
			next_invl()
			event = event - interval*(1. - noise)
			net_send(event, 1)
		}else if (w < 0 && on == 1) { 
			on = 0
		}
	}
	if (flag == 3) { 
		if (on == 0) {
			init_sequence(t)
			net_send(0, 1)
		}
	}
	if (flag == 1 && on == 1) {
		ispike = ispike + 1
		net_event(t)
		next_invl()
		if (on == 1) {
			net_send(event, 1)
		}
	}
}