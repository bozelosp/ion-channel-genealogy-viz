INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Gfluct
	RANGE g_e, E_e, g_e0,g_e1
	RANGE std_e, tau_e, D_e
	RANGE new_seed
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
	g_e0	= 0.0121 (umho)	
	std_e	= 0.0030 (umho)	
	tau_e	= 2.728	(ms)	
}

ASSIGNED {
	v	(mV)		
	i 	(nA)		
	g_e	(umho)		
	g_e1	(umho)		
	D_e	(umho umho /ms) 
	exp_e
	amp_e	(umho)
	xtemp	(1)
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
  i = g_e * (v - E_e)
}

PROCEDURE oup() {
  xtemp = normrand(0,1)
  if(tau_e==0) {
    g_e1 = std_e * xtemp
    g_e = g_e1
  } else {
    g_e1 = exp_e * g_e1 + amp_e * xtemp
    g_e = g_e0 + g_e1
  }
}

PROCEDURE new_seed(seed) {		
	set_seed(seed)
	VERBATIM
	  printf("Setting random generator with seed = %g\n", _lseed);
	ENDVERBATIM
}