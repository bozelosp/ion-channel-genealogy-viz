INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS orn
	RANGE g_e, g_e_max, cc_peak, g_e_baseline
	RANGE std_e, tau_e, D_e
	NONSPECIFIC_CURRENT i
        THREADSAFE 
        POINTER donotuse
}

UNITS {
	(nA) = (nanoamp) 
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	dt		  (ms)
        
	E_e	= 0 	  (mV)            
	g_e_max	= 75e-3 (umho)          
        cc_peak = 0     
        
        g_e_baseline    = 0       (umho)          
	std_e	= 1e-3    (umho)	  
	tau_e	= 400     (ms)            
}

ASSIGNED {
	v	(mV)		
	i 	(nA)		
	g_e	(umho)		
	g_e1	(umho)		
	D_e	(umho umho /ms) 
	exp_e
	amp_e	(umho)
        
        donotuse
}

STATE { O C D }

INITIAL {
	g_e1 = 0
	if(tau_e != 0) {
		D_e = 2 * std_e * std_e / tau_e
		exp_e = exp(-dt/tau_e)
		amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )
	}
        
        O = 0
        C = 1
        D = 0
}




BREAKPOINT {
        LOCAL SORN
	SOLVE oup
        SOLVE states METHOD derivimplicit
        
	if(tau_e==0) {
	   g_e = std_e * normrand123()
	}

        SORN = O * (1-D)
        
        g_e = g_e1 + SORN * cc_peak * g_e_max + g_e_baseline
          

        if(g_e < 0) {
            g_e = 0
        }

        i = g_e * (v - E_e)
}
DERIVATIVE states {
  LOCAL KO, KC1, KC2, KD1, KD2
  KO = 1/100
  KC1 = 1/100
  KC2 = 1e-4
  KD1 = 1/6000
  KD2 = 1/100
  O' = KO*(1-C-O)
  C' = KC1*(1-C)*C + KC2*(1-C)
  D' = KD1*O*(1-D) - KD2*D*(1-O)
}

PROCEDURE oup() {		
   if(tau_e!=0) {
	g_e1 =  exp_e * g_e1 + amp_e * normrand123()
   }
}

NET_RECEIVE(dummy) {
  C = 0
}


VERBATIM
double nrn_random_pick(void* r);
Rand* nrn_random_arg(int agpos);
ENDVERBATIM

FUNCTION normrand123() {
VERBATIM
	if (_p_donotuse) {
		
		_lnormrand123= nrn_random_pick(_p_donotuse);
	}else{
		
		if (_nt != nrn_threads) {
hoc_execerror("multithread random in NetStim"," only via hoc Random");
		}
ENDVERBATIM
		
		
		
		
		normrand123 = normrand(0,1)
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