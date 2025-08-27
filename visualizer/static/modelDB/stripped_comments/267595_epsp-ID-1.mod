INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	POINT_PROCESS epsp
	RANGE onset, tau0, tau1, imax, i, myv
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
	imax=0 	 (nA)
	v	 (mV)
}

ASSIGNED { i (nA)  myv (mV)}

LOCAL   a[2]
LOCAL   tpeak
LOCAL   adjust
LOCAL   amp

BREAKPOINT {
	myv = v
        i = curr(t)
}

FUNCTION myexp(x) {
	if (x < -100) {
	myexp = 0
	}else{
	myexp = exp(x)
	}
}

FUNCTION curr(x) {				
	tpeak=tau0*tau1*log(tau0/tau1)/(tau0-tau1)
	adjust=1/((1-myexp(-tpeak/tau0))-(1-myexp(-tpeak/tau1)))
	amp=adjust*imax
	if (x < onset) {
		curr = 0
	}else{
		a[0]=1-myexp(-(x-onset)/tau0)
		a[1]=1-myexp(-(x-onset)/tau1)
		curr = -amp*(a[0]-a[1])
	}
}