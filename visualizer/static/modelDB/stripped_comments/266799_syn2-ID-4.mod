INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS syn2
	RANGE onset, tau0, tau1, gmax, e, i, myv
	NONSPECIFIC_CURRENT i
}
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
}

PARAMETER {
	onset=0  (ms)
	tau0=0.2 (ms)
	tau1=3.0 (ms)
	gmax=0 	 (umho)
	e=0	 (mV)
	v	 (mV)
}

ASSIGNED { i (nA)  g (umho) myv (mV)}

LOCAL   a[2]
LOCAL   tpeak
LOCAL   adjust
LOCAL   amp

BREAKPOINT {
	myv = v
        g = cond(t)
	i = g*(v - e)
}

FUNCTION myexp(x) {
	if (x < -100) {
	myexp = 0
	}else{
	myexp = exp(x)
	}
}

FUNCTION cond(x) {
	tpeak=tau0*tau1*log(tau0/tau1)/(tau0-tau1)
	adjust=1/((1-myexp(-tpeak/tau0))-(1-myexp(-tpeak/tau1)))
	amp=adjust*gmax
	if (x < onset) {
		cond = 0
	}else{
		a[0]=1-myexp(-(x-onset)/tau0)
		a[1]=1-myexp(-(x-onset)/tau1)
		cond = amp*(a[0]-a[1])
	}
}