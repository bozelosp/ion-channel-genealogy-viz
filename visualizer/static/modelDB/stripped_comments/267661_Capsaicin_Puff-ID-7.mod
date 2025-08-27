NEURON {
	POINT_PROCESS Caps_Puff
	RANGE onset, tau_act, tau_inact, gmax,X, i, pump
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	onset=100 (ms)
	tau_act=1e6 (ms)	<1e-3,1e10>
	tau_inact=6500 (ms) <1e-3,1e10> 
	gmax=9 	(uS)	<0,1e9>
	X=1e-6 (M) <1e-15, 1>
	pump=500
	
	
	A1=-0.00717
	A2=1.00092
	x0=-7.32034
	dx=0.18614
	}

ASSIGNED { i (nA)  g (uS)}

BREAKPOINT {
	if (t<onset) {i=0}
	if (gmax) { at_time(onset) }
	if(t>=onset && t<onset+pump) {
	g = gmax *Conc(X)*alpha( (t - onset)/tau_act)
	i = -g
	}
	
	at_time(onset+pump) 
	if(t>=onset+pump) {
	g = gmax *Conc(X)*alpha( (pump)/tau_act) *beta((t - onset-pump)/tau_inact)
	i = -g
	}
}
FUNCTION alpha(x) {
	if (x < 0 || x > 10) {
		alpha = 0
	}else{
		alpha = (1-exp(-x))
	}
}

FUNCTION beta(x) {
	if (x < 0 || x > 10) {
		beta = 0
	}else{
		beta =exp(-x)
	}
}

FUNCTION Conc(x) {    
	
	
	Conc = A2 + (A1-A2)/(1 + exp((log10(x)-x0)/dx))
}