INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX Gfluctdv
	RANGE g_e, g_i, E_e, E_i, g_e0, g_i0, g_e1, g_i1
	RANGE std_e, std_i, tau_e, tau_i, D_e, D_i
	RANGE new_seed
	GLOBAL multex,multin
        NONSPECIFIC_CURRENT i

}

UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
   	(S) = (siemens)
}

PARAMETER {
	dt		(ms)

	E_e	= 0 	(mV)	
	E_i	= -75 	(mV)	


     	g_e0 = 0.0001 (S/cm2)	
     	g_i0 = 0.0005 (S/cm2)	

     	std_e = 3e-5 (S/cm2)	
     	std_i = 6e-5 (S/cm2)	

	tau_e	= 2.728	(ms)	
	tau_i	= 10.49	(ms)	

	multex=0
	multin=0
}

ASSIGNED {
	v	(mV)		
	i 	(mA/cm2)	
	g_e	(S/cm2)		
	g_i	(S/cm2)		
	g_e1	(S/cm2)		
	g_i1	(S/cm2)		
	D_e	(umho umho /ms) 
	D_i	(umho umho /ms) 
	exp_e
	exp_i
	amp_e	(umho)
	amp_i	(umho)
}

INITIAL {
	g_e1 = 0
	g_i1 = 0
	if(tau_e != 0) {
		D_e = 2 * std_e * std_e / tau_e
		exp_e = exp(-dt/tau_e)
		amp_e =sqrt(multex)*std_e * sqrt( (1-exp(-2*dt/tau_e)) )
	}
	if(tau_i != 0) {
		D_i = 2 * std_i * std_i / tau_i
		exp_i = exp(-dt/tau_i)
		amp_i = sqrt(multin)*std_i * sqrt( (1-exp(-2*dt/tau_i)) )
	}
}

BREAKPOINT {
	SOLVE oup
	if(tau_e==0) {
	   g_e = std_e * normrand(0,1)
	}
	if(tau_i==0) {
	   g_i = std_i * normrand(0,1)
	}
	g_e = multex*g_e0 + g_e1
	if(g_e < 0) { g_e = 0 }
	g_i = multin* g_i0 + g_i1
	if(g_i < 0) { g_i = 0 }
	i = g_e * (v - E_e) + g_i * (v - E_i)
}


PROCEDURE oup() {		
   if(tau_e!=0) {
	amp_e =sqrt(multex)*std_e * sqrt( (1-exp(-2*dt/tau_e)) )
	g_e1 =  exp_e * g_e1 + amp_e * normrand(0,1)
   }
   if(tau_i!=0) {
	amp_i = sqrt(multin)*std_i * sqrt( (1-exp(-2*dt/tau_i)) )
	g_i1 =  exp_i * g_i1 + amp_i * normrand(0,1)
   }
}


PROCEDURE new_seed(seed) {		
	set_seed(seed)
	VERBATIM
	  printf("Setting random generator with seed = %g\n", _lseed);
	ENDVERBATIM
}