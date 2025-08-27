NEURON {
	POINT_PROCESS DCNsynGABA
	NONSPECIFIC_CURRENT i
	RANGE g, i, e, tauRise, tauFall, startDeprLevel, deprLevel
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
}

PARAMETER {
	tauRise = 1 (ms)
	tauFall = 1 (ms)
	e = 0 (mV)
	startDeprLevel = 1 
	        
	        
}

ASSIGNED {
    relProbSS (1) 
    relProb[2] (1) 
    freq (1/s) 
    tau (ms)
    tSpikes[2] (ms)
    ISI (ms)
    deprLevel (1) 
                  
    notFirstSpike (1) 
                      

	v (mV)
	i (nA)
	g (microsiemens)
	factor
}

STATE {
	A (microsiemens)
	B (microsiemens)
}

INITIAL {
	LOCAL tp
	if (tauRise/tauFall > .9999) {
		tauRise = .9999*tauFall
	}
	A = 0
	B = 0
	tp = (tauRise*tauFall)/(tauFall - tauRise) * log(tauFall/tauRise)
	factor = -exp(-tp/tauRise) + exp(-tp/tauFall)
	factor = 1/factor

    notFirstSpike = 0
}

BREAKPOINT {
    
    
  	SOLVE state METHOD cnexp
   	g = (B - A) * deprLevel
   	i = g*(v - e)   	
}

DERIVATIVE state {
	A' = -A/tauRise
	B' = -B/tauFall
}

NET_RECEIVE(weight (microsiemens)) {
    deprLevel = giveFractionG()
    state_discontinuity(A, A + weight*factor)
	state_discontinuity(B, B + weight*factor)
}

FUNCTION giveFractionG() {
	if (notFirstSpike) {
        
        
        tSpikes[0] = tSpikes[1]
        tSpikes[1] = t
        ISI = tSpikes[1] - tSpikes[0]
        freq = 1000 / ISI
        
        relProbSS = 0.08 + 0.60*exp(-2.84*freq) + 0.32*exp(-0.02*freq)
        tau = 2 + 2500*exp(-0.274*freq) + 100*exp(-0.022*freq)
        relProb[1] = relProb[0] + (relProbSS - relProb[0]) * (1-exp(-ISI/tau))
        relProb[0] = relProb[1]

        giveFractionG = relProb[1]
	} else {
	    tSpikes[1] = t
	    relProb[0] = startDeprLevel
	    notFirstSpike = 1
	    giveFractionG = relProb[0]
	}
}