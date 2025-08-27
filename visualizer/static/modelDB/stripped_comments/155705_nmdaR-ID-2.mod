NEURON {
	POINT_PROCESS nmdaR
	RANGE g, g_eff, s, synon, t1, t2
	RANGE Cdur, Alpha, Beta, Erev, Rinf, Rtau
	RANGE mag, eta, gamma
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(mM) = (milli/liter)
	(uS) = (microsiemens)
}

PARAMETER {

	Cdur	= 4.0	(ms)		
	Alpha	= 0.3	(/ms)		
	Beta	= 0.01	(/ms)		
	Erev	= 0	    (mV)		
	mag     = 1     (mM)
	eta	    = 3.57  (mM)
	gamma   = 0.062 (/mV)
	g       = 1e-5	(uS)		
}


ASSIGNED {
	v		(mV)		
	i 		(nA)		
	g_eff 	(uS)		
	Rinf				
	Rtau	(ms)		
	synon
	t1
	t2
	s
}

STATE {Ron Roff}

INITIAL {
	Rinf = Alpha / (Alpha + Beta)
	Rtau = 1 / (Alpha + Beta)
	synon = 0
	t1 = 0
	t2 = 0
}

BREAKPOINT {
	SOLVE release METHOD cnexp
	s = Ron + Roff
	g_eff = g*s/(1 + mag * exp( - gamma * v ) / eta )
	i = g_eff * (v - Erev) 
}

DERIVATIVE release {
	Ron' = (synon*Rinf - Ron)/Rtau
	Roff' = -Beta*Roff
}







NET_RECEIVE(weight, on, nspike, r0, t0 (ms) ) {
	
        if (flag == 0) { 
		nspike = nspike + 1
		t1 = t1 + 1
		if (!on) {
			t2 = t2 + 1
			r0 = r0*exp(-Beta*(t - t0))
			t0 = t
			on = 1
			synon = synon + weight
			
			
			Ron=Ron + r0
			Roff=Roff - r0
		}
		
		net_send(Cdur, nspike)
        }
	if (flag == nspike) { 
		r0 = weight*Rinf + (r0 - weight*Rinf)*exp(-(t - t0)/Rtau)
		t0 = t
		synon = synon - weight
		
		
		Ron= Ron - r0
		Roff= Roff + r0
		on = 0
	}
}