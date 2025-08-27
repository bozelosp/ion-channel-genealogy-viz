NEURON {
	POINT_PROCESS AMPA_S
	RANGE g, gmax, synon
	NONSPECIFIC_CURRENT i
	GLOBAL Cmax, Cdur, Alpha, Beta, Erev, Rinf, Rtau
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	dt		(ms)
	deadtime = 1 (ms)	
	Cmax	= 0.5	(mM)		
	Cdur	= 0.3	(ms)		
	Alpha	= 0.94	(/ms mM)	
	Beta	= 0.18	(/ms)		
	Erev	= 0	(mV)		
	gmax		(umho)		
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g 		(umho)		
	Rinf				
	Rtau		(ms)		
	synon 			

}

STATE { Ron Roff }	
				
				

INITIAL {
	synon = 0
	Rinf = Alpha*Cmax / (Alpha*Cmax + Beta)
	Rtau = 1 / (Alpha*Cmax + Beta)
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


NET_RECEIVE(weight, on, r0, t0 (ms), tmp)  {
	if (flag == 0) {
		
		if (!on) {
			
			synon = synon + weight
			tmp = r0*exp(-Beta*(t-dt-t0)) 
			r0 = r0*exp(-Beta*(t-t0))
			Ron = Ron + tmp
			Roff = Roff - r0
			t0 = t 	
			on = 1
			net_send(Cdur,1)
		}
		
	}
	if (flag == 1) {
		
		synon = synon - weight
		
		tmp = weight*Rinf + (r0-weight*Rinf)*exp(-(t-dt-t0)/Rtau)	
		r0 = weight*Rinf + (r0-weight*Rinf)*exp(-(t-t0)/Rtau)
		Ron = Ron - r0
		Roff = Roff + tmp
		t0 = t		
		net_send(deadtime,2)	
	}
	if (flag == 2) {
		on = 0 
	}
}