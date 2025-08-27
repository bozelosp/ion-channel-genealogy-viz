NEURON	{ 
  ARTIFICIAL_CELL NetStimBox
  RANGE start, forcestop, status, nspk
  THREADSAFE 
  POINTER donotuse
}

PARAMETER {
	start		= 50 (ms)	
	forcestop 	= 200 (ms)	
	status 		= 0		
	nspk		= 1		
}

ASSIGNED {
	event (ms)
	on
	ispike
	donotuse
}


INITIAL {			
	on = 0  
}	

FUNCTION invl(mean (ms)) (ms) {				      
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

NET_RECEIVE (w) {
      
 
      if (status == 1) {
	
	FROM i=1 TO nspk {
  	      event = erand()*(forcestop-start)+start
        	
        	if (event < 0) {
        	         event = 0
        	}
		printf ("NetStimBox
		net_event(event)
	}
	status = 0			
      } 
 
}