INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS GABAa_S
	RANGE C, R, R0, R1, g, Cmax
	NONSPECIFIC_CURRENT i
	GLOBAL Cdur, Alpha, Beta, Erev, blockTime, Rinf, Rtau
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cmax	= 0.5	(mM)		
	Cdur	= 0.3	(ms)		
	Alpha	= 10.5	(/ms mM)	
	Beta	= 0.166	(/ms)		
	Erev	= -80	(mV)		
	blockTime = 2	(ms)		
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	C		(mM)		
	R				
	R0				
	Rinf				
	Rtau		(ms)		
	on				
	gmax				
	tLast
	nspike
	collisionBlock
}

INITIAL {
	R = 0
	C = 0
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / ((Alpha * Cmax) + Beta)
	on = 0
	R0 = 0
	nspike = 0
	collisionBlock = 0
}

BREAKPOINT {
	SOLVE release
	i = R*(v - Erev)
}

PROCEDURE release() {


	if (on) {				

	   R = gmax*Rinf + (R0 - gmax*Rinf) * exptable (- (t - tLast) / Rtau)
				
	} else {				

  	   R = R0 * exptable (- Beta * (t - tLast))
	}

	VERBATIM
	return 0;
	ENDVERBATIM
}








NET_RECEIVE(weight, ncType, ncPrb) {LOCAL ok, tmp

	
	

	INITIAL {
	}

	
      if (flag == 0) { 
		ok = 0

		if (ncType == 1) {
			collisionBlock = collisionBlock + 1
			net_send(blockTime, -1)





			ok = 1

		}
		else 
		if (collisionBlock == 0) {
			ok = 1
		}

		if (ok) {
			if (!on) {
				on = 1
				tLast = t
				R0 = R
				gmax = weight	
			}

			nspike = nspike + 1
			
			net_send(Cdur, nspike)			
		}
      }
	else
	if (flag == nspike) { 
		if (on) {
			on = 0
			tLast = t
			R0 = R
			gmax = 0
		}
	} 
	else
	if (flag == -1) {
		collisionBlock = collisionBlock - 1
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