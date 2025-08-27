NEURON	{ 
  ARTIFICIAL_CELL ThetaStim
  RANGE interval, number, start, actual_start
  RANGE noise
  RANGE outer_interval, outer_number, outer_start
  RANGE outer_noise
  THREADSAFE 
  POINTER donotuse
}

PARAMETER {
	interval	= 25 (ms) <1e-9,1e9>
	number	= 4 <0,1e9>	
	start		= 25 (ms)	
	noise		= 0 <0,1>	
	outer_interval	= 200 (ms) <1e-9,1e9>
	outer_number	= 5 <0,1e9>	
	outer_start		= 25 (ms)	
	outer_noise		= 0 <0,1>	
}

ASSIGNED {
	event (ms)
	on
	ispike
	donotuse
	outer_event (ms)
	outer_on
	outer_ispike
	actual_start 
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
	if (outer_noise < 0) {
		outer_noise = 0
	}
	if (outer_noise > 1) {
		outer_noise = 1
	}
	if ((outer_start >= 0 && outer_number > 0) && (number > 0)) {
		outer_on = 1
		
		
		outer_event = outer_start + invl(outer_interval) - outer_interval*(1. - outer_noise)
		
		if (outer_event < 0) {
			outer_event = 0
		}
		net_send(outer_event, 13) 
	}
}	

PROCEDURE init_sequence(t(ms)) {
	if (number > 0) {
		on = 1
		event = 0
		ispike = 0
	}
}
PROCEDURE init_outer_sequence(t(ms)) {
	if (outer_number > 0) {
		outer_on = 1
		outer_event = 0
		outer_ispike = 0
	}
}

FUNCTION invl(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) 
	}
	if (noise == 0) {
		invl = mean
	}else{
		invl = (1. - noise)*mean + noise*mean*erand()
	}
}
FUNCTION outer_invl(mean (ms)) (ms) {
	if (mean <= 0.) {
		mean = .01 (ms) 
	}
	if (outer_noise == 0) {
		outer_invl = mean
	}else{
		outer_invl = (1. - outer_noise)*mean + outer_noise*mean*erand()
	}
}
VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION erand() {
VERBATIM
	if (_p_donotuse) {
		
		_lerand = nrn_random_pick(_p_donotuse);
	}else{
		
		if (_nt != nrn_threads) {
hoc_execerror("multithread random in NetStim"," only via hoc Random");
		}
ENDVERBATIM
		
		
		
		
		erand = exprand(1)
VERBATIM
	}
ENDVERBATIM
}

PROCEDURE noiseFromRandom() {
VERBATIM
 {
	void** pv = (void**)(&_p_donotuse);
	if (ifarg(1)) {
		*pv = nrn_random_arg(1);
	}else{
		*pv = (void*)0;
	}
 }
ENDVERBATIM
}

PROCEDURE next_invl() {
	if (number > 0) {
		event = invl(interval)
	}
	if (ispike >= number) {
		on = 0
	}
}
PROCEDURE next_outer_invl() {
	if (outer_number > 0) {
		outer_event = outer_invl(outer_interval)
	}
	if (outer_ispike >= outer_number) {
		outer_on = 0
	}
}

NET_RECEIVE (w) {
	if (flag == 0) { 
		if (w > 0 && on == 0) { 
			
			init_sequence(t)
			
			
			next_invl()
			event = event - interval*(1. - noise)
			net_send(event, 1)
		}else if (w < 0) { 
			on = 0
		}
	}
	if (flag == 3) { 
		if (on == 1) { 
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
	if (flag == 13) { 
		if (outer_on == 1) { 
			init_outer_sequence(t)
			on = 1
			net_send(0, 11)           
		}
	}
	if (flag == 11 && outer_on == 1) {
		on = 1
		net_send(start, 3) 
		outer_ispike = outer_ispike + 1
		
		next_outer_invl()
		if (outer_on == 1) {
			net_send(outer_event, 11) 
		}
	}
}