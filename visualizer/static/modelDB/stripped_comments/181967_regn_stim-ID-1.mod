NEURON	{ 
  ARTIFICIAL_CELL RegnStim
  RANGE interval, number, start
  RANGE noise
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
	tspike	
	donotuse
}

PROCEDURE seed(x) {
	set_seed(x)
}

INITIAL {
	on = 0 
	tspike = start
	ispike = 0
	if (noise < 0) {
		noise = 0
	}
	if (noise > 1) {
		noise = 1
	}
	if (start >= 0 && number > 0) {
		on = 1
		
		event = start + noise*interval*erand()
		
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

		invl = tspike + mean + noise*mean*erand() - t
		if (invl <= 0) {
			invl = .01 (ms)	
		}



	}
	tspike = tspike + mean
}
VERBATIM
#ifndef NRN_VERSION_GTEQ_8_2_0
double nrn_random_pick(void* r);
void* nrn_random_arg(int argpos);
#define RANDCAST
#else
#define RANDCAST (Rand*)
#endif

ENDVERBATIM

FUNCTION erand() {
VERBATIM
	if (_p_donotuse) {
		
		_lerand = nrn_random_pick(RANDCAST _p_donotuse);
	}else{
ENDVERBATIM
		
		
		
		
		erand = normrand(0, 1)
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