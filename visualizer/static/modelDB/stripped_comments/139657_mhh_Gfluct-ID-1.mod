INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS mhh_Gfluct
	RANGE g_e, E_e, g_e0, g_e1
	RANGE std_e, tau_e, D_e
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp) 
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	dt		(ms)

	E_e	= 0 	(mV)	
	g_e0	= 0.000001 (umho)	
	std_e	= 0.0002 (umho)	
	tau_e	= 2	(ms)	
}

ASSIGNED {
	v	(mV)		
	i 	(nA)		
	g_e	(umho)		
	g_e1	(umho)		
	D_e	(umho umho /ms) 
	exp_e
	amp_e	(umho)
}

INITIAL {
	g_e1 = 0
	if(tau_e != 0) {
		D_e = 2 * std_e * std_e / tau_e
		exp_e = exp(-dt/tau_e)
		amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )
	}
}

BREAKPOINT {
	SOLVE oup
	if(tau_e==0) {
	   g_e = std_e * normrand(0,1)
	}
	g_e = g_e0 + g_e1
	i = g_e * (v - E_e)
}


PROCEDURE oup() {		
   if(tau_e!=0) {
	g_e1 =  exp_e * g_e1 + amp_e * normrand(0,1)
   }
}


PROCEDURE new_seed(seed) {		
	set_seed(seed)
	VERBATIM
	  printf("Setting random generator with seed = %g\n", _lseed);
	ENDVERBATIM
}