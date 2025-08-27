NEURON	{ 
  ARTIFICIAL_CELL pNetStim 
  RANGE interval, number, start
  RANGE noise
  THREADSAFE 
  POINTER donotuse
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
	donotuse
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
        printf("INITIAL t=%g\n", t)
	on = 0 
	ispike = 0
	if (noise < 0) {
		noise = 0
	}
	if (noise > 1) {
		noise = 1
	}
	if (start >= 0 && number > 0) {
		on = 1
		
		
		event = start + invl(interval) - interval*(1. - noise)
		
		if (event < 0) {
			event = 0
		}
		net_send(event, 3)
		printf("initializing net_send with event = %g\n", event)
	}
}	

PROCEDURE init_sequence(t(ms)) {
	printf("init_sequence with number = %g, ignored arg t = %g\n", number, t)
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
		invl = (1. - noise)*mean + noise*mean*erand()
	}
	printf("invl set to  = %g at t=%g\n", invl, t)
}
VERBATIM
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
ENDVERBATIM

FUNCTION erand() {
VERBATIM
	if (_p_donotuse) {
		
		_lerand = nrn_random_pick(_p_donotuse);
		printf("_lerand = %g\n", _lerand);
	}else{
		
		if (_nt != nrn_threads) {
hoc_execerror("multithread random in NetStim"," only via hoc Random");
		}
ENDVERBATIM
		
		
		
		
		erand = exprand(1)
		printf("erand = %g\n", erand)
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
	printf("next_invl
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
	
}