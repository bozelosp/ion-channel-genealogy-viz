NEURON {
	POINT_PROCESS AMPA
	RANGE R, gmax, g, ina, Alpha, Beta, iAMPA
	USEION na WRITE ina
	NONSPECIFIC_CURRENT  iAMPA
	GLOBAL Cdur, Erev, Rinf, Rtau
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
    Cmax	= 0.1	(mM)		

	Cdur	= 1.1	(ms)		

	Alpha	= 1	(/ms)	

	Beta	= 0.5 (/ms)		
	Erev	= 0	(mV)		
	gmax    = 1  (uS)
}


ASSIGNED {
	v		(mV)		
	iAMPA 		(nA)		
	g 		(uS)		
	Rinf				
	Rtau		(ms)		
	synon
	ina
	ica
}

STATE {Ron Roff}

INITIAL {
    Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
    Rtau = 1 / ((Alpha * Cmax) + Beta)
	synon = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
	g = (Ron + Roff)* gmax
	iAMPA = g*(v - Erev)
	ina = 0.9*iAMPA
	iAMPA = 0.1*iAMPA

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