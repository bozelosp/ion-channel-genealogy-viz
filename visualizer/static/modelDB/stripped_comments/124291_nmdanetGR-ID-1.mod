NEURON {
	POINT_PROCESS nmdanet
	RANGE R, g, mg, Alpha, Beta, nmdafactor
	NONSPECIFIC_CURRENT i
	GLOBAL Cdur, Erev, Rinf, Rtau
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cdur	= 1		(ms)	
	Alpha	= 0.1		(/ms)	
	Beta	= 0.5		(/ms)	


	Erev	= 0	(mV)		
	mg	= 1    (mM)		
	nmdafactor=2.65
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	Rinf				
	Rtau		(ms)		
	synon
}

STATE {Ron Roff}

INITIAL {
	Rinf = Alpha / (Alpha + Beta)
	Rtau = 1 / (Alpha + Beta)
	synon = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
	g = nmdafactor*mgblock(v)*(Ron + Roff)*1(umho)
	i = g*(v - Erev)		
}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
}








FUNCTION mgblock(v(mV)) {
	TABLE 
	DEPEND mg
	FROM -140 TO 80 WITH 1000

	

	mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}


NET_RECEIVE(weight, on, nspike, r0, t0 (ms)) {
	
        if (flag == 0) { 
		nspike = nspike + 1
		if (!on) {
			r0 = r0*exp(-Beta*(t - t0))
			t0 = t
			on = 1
			synon = synon + weight
			state_discontinuity(Ron, Ron + r0)
			state_discontinuity(Roff, Roff - r0)
		}
		
		net_send(Cdur, nspike)
        }
	if (flag == nspike) { 
		r0 = weight*Rinf + (r0 - weight*Rinf)*exp(-(t - t0)/Rtau)
		t0 = t
		synon = synon - weight
		state_discontinuity(Ron, Ron - r0)
		state_discontinuity(Roff, Roff + r0)
		on = 0
	}
}