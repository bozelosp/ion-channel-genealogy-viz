VERBATIM


static char rcsid[] = "$Id: rglu_score.mod,v 1.4 2000/08/14 22:21:27 karchie Exp $";

extern double exprand(double mean);

ENDVERBATIM


INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS rGluSc
	RANGE C, gp_NMDA, lastrelease, score, idx, phase, mean_ia_time, gmax_NMDA, R_NMDA, R0_NMDA, R1_NMDA, gp_AMPA, gmax_AMPA, R_AMPA, R0_AMPA, R1_AMPA
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur, Deadtime, score_thresh, score_tau, Alpha_NMDA, Beta_NMDA, Erev_NMDA, Rinf_NMDA, Rtau_NMDA, Alpha_AMPA, Beta_AMPA, Erev_AMPA, Rinf_AMPA, Rtau_AMPA
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt (ms)
	
	Cmax	= 1	(mM)		
	Cdur	= 1.1	(ms)		
	Deadtime = 0	(ms)		
	phase	= 0			
	mean_ia_time = 100 (ms)		

	
	Alpha_NMDA = 10	(/ms mM)	
	Beta_NMDA = 0.0125 (/ms)	
	Alpha_AMPA = 10 (/ms mM)	
	Beta_AMPA = 0.5 (/ms)		

	
	Erev_NMDA	= 0	(mV)		
	Erev_AMPA = 0	(mV)

	
	gmax_NMDA	(umho)		
	gmax_AMPA	(umho)		

	
	score_thresh = -63	(mV)	
        score_tau = 20	(ms)		

	
	eta     = 0.33  (/mM)
	mag     = 1     (mM)
	gamma   = 0.06  (/mV)
}

ASSIGNED {
	v		(mV)		
	i 		(nA)		
	gp_NMDA		(umho)		
	gp_AMPA		(umho)		

	depol		(mV)		
	C		(mM)		

	R_NMDA				
	R0_NMDA				
	R1_NMDA				
	R_AMPA				
	R0_AMPA				
	R1_AMPA				

	Rinf_NMDA			
	Rtau_NMDA	(ms)		
	Rinf_AMPA			
	Rtau_AMPA	(ms)		

	lastrelease	(ms)		
	ia_time		(ms)		
	score				
	score_postsyn_act		
	idx				
}

INITIAL {
	C = 0

	R_NMDA = 0
	R1_NMDA = R_NMDA
	R_AMPA = 0
	R1_AMPA = R_AMPA

	Rinf_NMDA = Cmax*Alpha_NMDA / (Cmax*Alpha_NMDA + Beta_NMDA)
	Rtau_NMDA = 1 / ((Alpha_NMDA * Cmax) + Beta_NMDA)

	Rinf_AMPA = Cmax*Alpha_AMPA / (Cmax*Alpha_AMPA + Beta_AMPA)
	Rtau_AMPA = 1 / ((Alpha_AMPA * Cmax) + Beta_AMPA)

	score = 0
	score_postsyn_act = 0

	
	
	
VERBATIM
	assert(phase < 1.0);
	if (phase < 0.0) {
		assert(mean_ia_time > 0.0);
	} else {
		assert(mean_ia_time > Cdur + Deadtime);
	}
ENDVERBATIM

	
	if (phase < 0.0) {
		
		lastrelease = t - next_ia_time()
		ia_time = next_ia_time()
	} else {
		
		ia_time = mean_ia_time
		lastrelease = t - phase * mean_ia_time
	}	
}

BREAKPOINT {
	SOLVE release
	gp_AMPA = gmax_AMPA * R_AMPA
	gp_NMDA = (gmax_NMDA * R_NMDA)/(1 + eta * mag * exp( - (gamma * v)))
	i = gp_AMPA*(v - Erev_AMPA) + gp_NMDA*(v - Erev_NMDA)

	
	depol = v - score_thresh	
	if (depol < 0) {
		depol = 0	
	}

	
	score_postsyn_act = score_postsyn_act + dt * (depol - score_postsyn_act) / score_tau

	score = score + dt * gp_NMDA * score_postsyn_act
}

PROCEDURE release() { LOCAL q
	
	
	q = t - lastrelease
	if (q >= ia_time) {			
		C = Cmax
		R0_NMDA = R_NMDA
		R0_AMPA = R_AMPA
		lastrelease = t
		ia_time = next_ia_time()	
	} else if (q < Cdur) {			
		
	} else if (C == Cmax) {			
		R1_NMDA = R_NMDA
		R1_AMPA = R_AMPA
		C = 0.
	}		
	
	if (C > 0) {				
	   R_NMDA = Rinf_NMDA + (R0_NMDA - Rinf_NMDA) * exptable (- (t - lastrelease) / Rtau_NMDA)
	   R_AMPA = Rinf_AMPA + (R0_AMPA - Rinf_AMPA) * exptable (- (t - lastrelease) / Rtau_AMPA)
	} else {				
  	   R_NMDA = R1_NMDA * exptable(-Beta_NMDA * (t - (lastrelease + Cdur)))
  	   R_AMPA = R1_AMPA * exptable(-Beta_AMPA * (t - (lastrelease + Cdur)))
	}

	VERBATIM
	return 0;
	ENDVERBATIM
}


FUNCTION next_ia_time() {
	if (phase < 0.0) {
		
		next_ia_time = 0
		while (next_ia_time <= Cdur + Deadtime) {
VERBATIM
			_lnext_ia_time =
				exprand(mean_ia_time);
ENDVERBATIM
		}
	} else {
		
		next_ia_time = mean_ia_time
	}
}


FUNCTION exptable(x) { 
	TABLE  FROM -10 TO 10 WITH 2000

	if ((x > -10) && (x < 10)) {
		exptable = exp(x)
	} else {
		exptable = 0.
	}
}