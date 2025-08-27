NEURON {
	POINT_PROCESS NMDA
	RANGE g, Alpha, Beta, e, gmax, ica, Cdur, iNMDA
	USEION ca WRITE ica
	NONSPECIFIC_CURRENT  iNMDA
	GLOBAL mg, Cmax
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
	(mM) = (milli/liter)
}

PARAMETER {
	Cmax	= 1	 (mM)           
	Cdur	= 1	 (ms)		
	Alpha	= 4 (/ms /mM)	
	Beta 	=0.01   (/ms)
	e	= 0	 (mV)		
    mg   = 1      (mM)           
	gmax = 1   (uS)

}


ASSIGNED {
	v		(mV)		
	iNMDA 		(nA)		
	g 		(uS)		
	Rinf				
	Rtau		(ms)		
	synon
    B                       
	ica
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
	g = (Ron + Roff)* gmax * B
	iNMDA = g*(v - e)
    ica = 7*iNMDA/10   
    iNMDA = 3*iNMDA/10

}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
}

FUNCTION mgblock(v(mV)) {
        TABLE
        DEPEND mg
        FROM -140 TO 80 WITH 1000

        


	 mgblock = 1 / (1 + exp(0.072 (/mV) * -v) * (mg / 3.57 (mM)))  

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