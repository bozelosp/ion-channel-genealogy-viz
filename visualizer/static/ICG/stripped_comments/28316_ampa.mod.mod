NEURON {
	POINT_PROCESS AMPA_S
	NONSPECIFIC_CURRENT i
	RANGE R, g, gmax, i
	GLOBAL Cdur_a, Alpha_a, Beta_a, Erev_a, Rinf_a, Rtau_a
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {

	Cdur_a	= 1	(ms)		
	Alpha_a	= 1.1	(/ms)	
	Beta_a	= 0.19	(/ms)		
	Erev_a	= 0	(mV)		
	gmax 
}


ASSIGNED {
	v		(mV)		
	i		(nA)		
	g 		(umho)		
	Rinf_a				
	Rtau_a		(ms)		
	synon
}

STATE {Ron Roff}

INITIAL {
	Rinf_a = Alpha_a / (Alpha_a + Beta_a)
	Rtau_a = 1 / (Alpha_a + Beta_a)
	synon = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
	g = gmax*(Ron + Roff)*1(umho)
	i = g*(v - Erev_a)
}

DERIVATIVE release {
	Ron' = (synon*Rinf_a - Ron)/Rtau_a
	Roff' = -Beta_a*Roff
}







NET_RECEIVE(weight, on, nspike, r0, t0 (ms)) {
	
        if (flag == 0) { 
		nspike = nspike + 1
		if (!on) {
			r0 = r0*exp(-Beta_a*(t - t0))
			t0 = t
			on = 1
			synon = synon + weight
			state_discontinuity(Ron, Ron + r0)
			state_discontinuity(Roff, Roff - r0)
		}
		
		net_send(Cdur_a, nspike)
        }
	if (flag == nspike) { 
		r0 = weight*Rinf_a + (r0 - weight*Rinf_a)*exp(-(t - t0)/Rtau_a)
		t0 = t
		synon = synon - weight
		state_discontinuity(Ron, Ron - r0)
		state_discontinuity(Roff, Roff + r0)
		on = 0
	}
}