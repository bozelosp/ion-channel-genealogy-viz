INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX NMDAKIT
	NONSPECIFIC_CURRENT i
	USEION ca READ cai WRITE ica
	RANGE onset,period, nbre, tau0, tau1, e, B, cao, gmax, g
        GLOBAL Erev, mg, temp, F, R
}
UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
        (mM) = (milli/liter)
        (celsius) = (degC)
	
}

PARAMETER {
	onset = 10  (ms)
	period = 50 (ms)
	nbre=20
	tau0 = 2.0 (ms)
	tau1 = 26.0 (ms)

	cao = 1.5	(mM)		
	cai		(mM)		
        Erev    = 0     (mV)            
        gmax           (S/cm2)          
        mg      = 1    (mM)             
	Px=4.6925	(cm3 mV/coulomb)	
	F = 96.49 (kilocoulomb)
	R = 8.314 (joule/degC)
	temp = 37		(degC)
}

ASSIGNED {
	ica             (mA/cm2)        
        v               (mV)            
        i               (mA/cm2)        
        g               (S/cm2)          
        B                               
}

LOCAL   a[2]
LOCAL   tpeak
LOCAL   adjust
LOCAL   amp

BREAKPOINT {
        B = mgblock(v)          
	g = cond(t,onset)
	if (nbre>1) {
	  FROM j=1 TO (nbre-1) {
	    g=g+cond(t,onset+j*period)
	  }
	}
	g=g*B
        ica = (0.001) * g * (0.051(cm3/coulomb)*v+Px)   *   (4.0*v*F*F / (R * (temp+273) ))   *   (-cao*exptable(-2*v*F/(R* (temp+273) ))  +  cai)  /  (1.0 - exptable(-2.0*v*F/(R* (temp+273) ))) 
	i = g*(v - Erev) - ica
}

FUNCTION myexp(x) {
	if (x < -100) {
	myexp = 0
	}else{
	myexp = exp(x)
	}
}

FUNCTION cond(x (ms), onset1 (ms)) (S/cm2) {
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


FUNCTION exptable(x) { 
        TABLE  FROM -10 TO 10 WITH 2000

        if ((x > -10) && (x < 10)) {
                exptable = exp(x)
        } else {
                exptable = 0.
        }
}

FUNCTION mgblock(v(mV)) {
        TABLE 
        DEPEND mg
        FROM -140 TO 80 WITH 1000

        

        mgblock = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}