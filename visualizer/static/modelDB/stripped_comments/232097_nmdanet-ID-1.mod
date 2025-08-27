NEURON {
	POINT_PROCESS NMDA
	RANGE g, Alpha, Beta, e
	NONSPECIFIC_CURRENT i
	GLOBAL Cdur, mg, Cmax
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	Cmax	= 1	 (mM)           
	Cdur	= 30	 (ms)		
	Alpha	= 0.072	 (/ms /mM)	
	Beta	= 0.0066 (/ms)		
	e	= 45	 (mV)		
        mg      = 1      (mM)           
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	Rinf				
	Rtau		(ms)		
	synon
        B                               
}

STATE {Ron Roff}

INITIAL {
	Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
	Rtau = 1 / (Cmax*Alpha + Beta)
	synon = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
        B = mgblock(v)
	g = (Ron + Roff)*1(umho) * B
	i = g*(v - e)
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