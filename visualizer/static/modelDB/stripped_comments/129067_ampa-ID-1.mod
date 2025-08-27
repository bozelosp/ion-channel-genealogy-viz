NEURON {
	POINT_PROCESS AMPA_S
	RANGE R, g
	NONSPECIFIC_CURRENT i
	GLOBAL Cdur, Alpha, Beta, Erev, Rinf, Rtau
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cdur	= 0.3	(ms)		
	Alpha	= 0.94	(/ms)	
	Beta	= 0.18	(/ms)		
	Erev	= 0	(mV)		
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
	g = (Ron + Roff)
	i = g*(v - Erev)
}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
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