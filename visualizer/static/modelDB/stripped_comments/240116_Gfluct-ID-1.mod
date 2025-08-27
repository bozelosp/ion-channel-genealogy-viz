INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS Gfluct
	RANGE g_e, g_i, E_e, E_i, g_e0, g_i0, g_e1, g_i1
	RANGE std_e, std_i, tau_e, tau_i, D_e, D_i
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
	dt		(ms)

	E_e	= 0 	(mV)	
	E_i	= -75 	(mV)	

	g_e0	= 0.0121 (umho)	
	g_i0	= 0.0573 (umho)	

	std_e	= 0.0030 (umho)	
	std_i	= 0.0066 (umho)	

	tau_e	= 3	(ms)	
	tau_i	= 10.49	(ms)	
}

ASSIGNED {
	v	(mV)		
	i 	(nA)		
	g_e	(umho)		
	g_i	(umho)		
	g_e1	(umho)		
	g_i1	(umho)		
	D_e	(umho umho /ms) 
	D_i	(umho umho /ms) 
	exp_e
	exp_i
	amp_e	(umho)
	amp_i	(umho)

        donotuse
}

INITIAL {
	g_e1 = 0
	g_i1 = 0
	if(tau_e != 0) {
		D_e = 2 * std_e * std_e / tau_e
		exp_e = exp(-dt/tau_e)
		amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )
	}
	if(tau_i != 0) {
		D_i = 2 * std_i * std_i / tau_i
		exp_i = exp(-dt/tau_i)
		amp_i = std_i * sqrt( (1-exp(-2*dt/tau_i)) )
	}
}

BREAKPOINT {
	SOLVE oup
	if(tau_e==0) {
	   g_e = std_e * normrand123()
	}
	if(tau_i==0) {
	   g_i = std_i * normrand123()
	}
	g_e = g_e0 + g_e1
	if(g_e < 0) { g_e = 0 }
	g_i = g_i0 + g_i1
	if(g_i < 0) { g_i = 0 }
	i = g_e * (v - E_e) + g_i * (v - E_i)
}


PROCEDURE oup() {		
   if(tau_e!=0) {
	g_e1 =  exp_e * g_e1 + amp_e * normrand123()
   }
   if(tau_i!=0) {
	g_i1 =  exp_i * g_i1 + amp_i * normrand123()
   }
}


PROCEDURE new_seed(seed) {		
	set_seed(seed)

}


VERBATIM
double nrn_random_pick(void* r);
Rand* nrn_random_arg(int argpos);
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