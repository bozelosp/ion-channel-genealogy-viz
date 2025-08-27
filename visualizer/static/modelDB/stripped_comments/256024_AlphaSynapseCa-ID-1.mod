NEURON {
 	POINT_PROCESS AlphaSynapseCa
 	RANGE onset, tau, gmax, e, i
    USEION ca READ ica, eca WRITE ica
 	NONSPECIFIC_CURRENT i
 }
 UNITS {
 	(nA) = (nanoamp)
 	(mV) = (millivolt)
 	(uS) = (microsiemens)
 }
 
 PARAMETER {
 	onset=0 (ms)
 	tau=0.1 (ms)	<1e-3,1e6>
 	gmax=0 	(uS)	<0,1e9>
 	cap=0.05 (1)
 	e=0	(mV)
 }
 
 ASSIGNED {
 	v (mV)
 	i (nA)
 	g (uS)
 	ica (nA)
 	eca (mV)
 }
 
 BREAKPOINT {
 	if (gmax) { at_time(onset) }
 	g = gmax * alpha( (t - onset)/tau )
 	i = g*(v - e)
 	ica = i*cap
 }
 
 FUNCTION alpha(x) {
 	if (x < 0 || x > 10) {
 		alpha = 0
 	}else{
 		alpha = x * exp(1 - x)
 	}
 }