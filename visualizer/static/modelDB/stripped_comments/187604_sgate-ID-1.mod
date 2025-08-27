NEURON  { 
  ARTIFICIAL_CELL SGate 
  RANGE period, number, start, phase
  RANGE depth, gid, randi
  THREADSAFE 
  POINTER donotuse
}

UNITS {
  PI = (pi) (1)
}

PARAMETER {
  period = 100 (ms) <1e-9,1e9>
  number = 1 <0,1e9> 
  start = 50 (ms) 
  depth = 0 <0,1> 
  phase = 0 (ms)
	gid = 0
	randi = 0
}

ASSIGNED {
  on (1)
  donotuse
  numtogo (1) 
  r (1)
}

INITIAL {
  if (period < 0) { period = 1e9 }
  if (number < 0) { number = 0 }
  if (start < 0) { start = 0 }
  if (phase < 0) { phase = 0 }
  if (depth < 0) { depth = 0 }
  if (depth > 1) { depth = 1 }

  on = 0 
  if (number > 0) {
    numtogo = number
    net_send(start, 1) 
  }
}  

PROCEDURE seed(x) {
  set_seed(x)
}

VERBATIM
#ifndef NRN_VERSION_GTEQ_8_2_0
double nrn_random_pick(void*);
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
    
    
    
    
    erand = scop_random()
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



FUNCTION p(t (ms)) {
  p = 0
  if (on == 1) {
    p = 1 + depth*(cos(2*PI*(t-phase)/period) - 1)/2
  }
}




NET_RECEIVE (w) {
  if (flag == 0) { 
    if (on == 1) {
      
        r = erand()
        if (r < p(t)) { net_event(t) }
    }
  } else if (flag == 1) {
    if (numtogo>0) { 
      on = 1
      numtogo = numtogo-1
      net_send(period, 1) 
    } else { 
      on = 0
    }
  }
}