COMMENT
Capsaicin point_process puff (Barkai et al. 2020, Goldstein et al. 2019)
ENDCOMMENT

COMMENT
an synaptic current with alpha function conductance defined by
        i = g * (v - e)      i(nanoamps), g(microsiemens);
        where
         g = 0 for t < onset and
         g = gmax * (t - onset)/tau * exp(-(t - onset - tau)/tau)
          for t > onset
this has the property that the maximum value is gmax and occurs at
 t = delay + tau.
ENDCOMMENT
			

COMMENT
A Puff pipette Current for Capsaicin 
        i = g * (v - e)      i(nanoamps), g(microsiemens);
        where
         g = 0 for t < onset and
         g = gmax *Con(X)* alpha*beta* for t > onset
		 
		 Alpha, Beta are the current open and close kinetics
		 Conc(X) is the logarthimic dose response curve 'gate', where X is the real concentration (not in logarithmic values)
ENDCOMMENT
			
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
	tau_inact=6500 (ms) <1e-3,1e10> :500
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

FUNCTION Conc(x) {    :Concentration 'gate'
	:Capsaicin dose-response curve, Boltzmann function fit according to Senning E,Gordon S et al. eLife 2015.
	:Note that the dose response curve accepts the concentration in Molars and 
	Conc = A2 + (A1-A2)/(1 + exp((log10(x)-x0)/dx))
}

