COMMENT
synaptic current with exponential rise and decay conductance defined by
        i = g * (v - e)      i(nanoamps), g(micromhos);
        where
         g = 0 for t < onset and
         g=amp*((1-exp(-(t-onset)/tau0))-(1-exp(-(t-onset)/tau1)))
          for t > onset 
ENDCOMMENT 
					       
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
     THREADSAFE
	POINT_PROCESS nmda
	RANGE onset, tau0, tau1, gmax, e, i, g
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	onset=2  (ms)
	tau0=3 (ms)
	tau1=90 (ms)
	gmax=.001	 (umho) 
	e=5	 (mV)
	v	 (mV)
	nmg = 0.3 (1) : eta*[Mg]o for eta = 0.25 (/mM) and [Mg]o = 1 mM
	  : 0.33 used in Polsky's glutamate.mod
	  : Zador et al. 1990 assumed eta = 0.33 (/mM)
	gamma = 0.08 (/mV) : .08 used in Polsky's glutamate.mod
	  : Zador et al. 1990 used 0.06 (/mV)
}

ASSIGNED { i (nA)  g (umho) }

LOCAL   a[2]
LOCAL   tpeak
LOCAL   adjust
LOCAL   amp

BREAKPOINT {
    g = cond(t)
	i = g*(v - e)
}

FUNCTION cond(x(ms))(umho) {
	tpeak=tau0*tau1*log(tau0/tau1)/(tau0-tau1)
	adjust=1/((1-exp(-tpeak/tau0))-(1-exp(-tpeak/tau1)))
	amp=adjust*gmax
	if (x < onset) {
		cond = 0
	}else{
		a[0]=1-exp(-(x-onset)/tau0)
		a[1]=1-exp(-(x-onset)/tau1)
		:cond = amp*(a[0]-a[1])
		cond = (amp*(a[0]-a[1])/(1+nmg*exp(-gamma*(v))))
	}
}
