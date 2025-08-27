NEURON {
	POINT_PROCESS FastInhib
	RANGE tau1, tau2, e, i
	NONSPECIFIC_CURRENT i
	RANGE gmax
	RANGE x, mgid, ggid, srcgid
	GLOBAL ltdinvl, ltpinvl, sighalf, sigslope

	RANGE g
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	tau1= 1 (ms) <1e-9,1e9> 
	tau2 = 200 (ms) <1e-9,1e9> 
	gmax = .003 (uS) 
	e = -80	(mV)
	ltdinvl = 250 (ms)		
	ltpinvl = 33.33 (ms)		
	sighalf = 25 (1)
	sigslope = 3 (1)
	x = 0 (um) 
	mgid = -1 
	ggid = -1 
	srcgid = -1 
}

ASSIGNED {
	v (mV)
	i (nA)
	g (uS)
	factor
	w (uS)
	total (uS)
}

STATE {
	A
	B
}

INITIAL {
	LOCAL tp
	if (tau1/tau2 > .9999) {
		tau1 = .9999*tau2
	}
	A = 0
	B = 0
	tp = (tau1*tau2)/(tau2 - tau1) * log(tau2/tau1)
	factor = -exp(-tp/tau1) + exp(-tp/tau2)
	factor = 1/factor
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	g = (B - A)*gmax
	i = g*(v - e)
}

DERIVATIVE state {
	A' = -A/tau1
	B' = -B/tau2
}

FUNCTION plast(step(1))(1) {
	plast = 1 - 1/(1 + exp((step - sighalf)/sigslope))
}

FUNCTION norm_weight_to_sig(w) {
    norm_weight_to_sig = floor(0.4999 + log(((-1/(w-1))-1)/exp(-sighalf/sigslope))*sigslope)
}

NET_RECEIVE(weight, s, w, tlast (ms)) {
	INITIAL {
		s = 0
		w = 0
		tlast = -1e9(ms)
	}
	if (t - tlast < ltpinvl) { 
		s = s + 1
		if (s > 2*sighalf) { s = 2*sighalf }
	}else if (t - tlast > ltdinvl) { 
	}else{ 
		s = s - 1
		if (s < 0) { s = 0 }
	}
	tlast = t
	w = weight 
	A = A + w*factor
	B = B + w*factor
}