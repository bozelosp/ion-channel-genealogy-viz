INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
        POINT_PROCESS AMPAKIT
	RANGE onset,periodpre, periodpost, delta,nbrepre, nbrepost, change, tau0, tau1, g,gmax, e, i,C
        NONSPECIFIC_CURRENT i
        GLOBAL Erev,Cmax,Cdur
}
UNITS {
        (celsius) = (degC)        
	(nA) = (nanoamp)
        (mV) = (millivolt)
        (umho) = (micromho)
        (mM) = (milli/liter)
}

PARAMETER {
	onset = 10  (ms)
	periodpre = 50 (ms)	
	periodpost=20
	delta= 10 (ms)		
	nbrepre=2			
	nbrepost=1
	tau0 = 0.34 (ms)
	tau1 = 2.0  (ms)
        Erev = 0    (mV)            
	Cmax	= 1	(mM)		
	Cdur	= 1	(ms)		
        gmax = 0.002 (umho)          

}


ASSIGNED {
        v               (mV)            
        i               (nA)            
        g               (umho)          
        C		(mM)		
	change
}

LOCAL   a[2]
LOCAL   tpeak
LOCAL   adjust
LOCAL   amp

INITIAL {
	C = 0
}


BREAKPOINT {
        g = cond(t,onset)
	C = trans(t,onset)
	if (nbrepre>1) {
	  FROM j=1 TO (nbrepre-1) {
	    g = g+cond(t,onset+j*periodpre)
	    C = C+trans(t,onset+j*periodpre)
	  }
	}
        i = g*(v - Erev)
}

FUNCTION myexp(x) {
	if (x < -100) {
	myexp = 0
	}else{
	myexp = exp(x)
	}
}

FUNCTION cond(x (ms), onset1 (ms)) (umho) {
	tpeak=tau0*tau1*log(tau0/tau1)/(tau0-tau1)
	adjust=1/((1-myexp(-tpeak/tau0))-(1-myexp(-tpeak/tau1)))
	amp=adjust*gmax
	if (x < onset1) {
		cond = 0
	}else{
		a[0]=1-myexp(-(x-onset1)/tau0)
		a[1]=1-myexp(-(x-onset1)/tau1)
		cond = amp*(a[0]-a[1])
	}
}
FUNCTION trans(x (ms), onset1 (ms)) (mM) {
	if ((x>onset1) && (x-onset1<=Cdur)) {
		trans=Cmax
	} else {
		trans=0
	}
}