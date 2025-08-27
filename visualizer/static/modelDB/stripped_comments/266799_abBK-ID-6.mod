NEURON {
	SUFFIX abBK
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE gabkbar,gabk, ik

}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)

}
CONSTANT {
    q10 = 3
}

PARAMETER {
	gabkbar = .01	(S/cm2)	
    cai (mM)
    base = 1  	(mV)	
}

ASSIGNED {
	v		(mV)
	ek		(mV)
	ik		(mA/cm2)
    gabk		(S/cm2)
    abinf		(mV)
    abtau		(ms)
    qt
}

STATE { ab }

BREAKPOINT {
	SOLVE state METHOD cnexp
	gabk = gabkbar*ab
	ik = (gabk)*(v - ek)
}

DERIVATIVE state {	
	rates(v, cai)				      
	ab' = (abinf-ab)/abtau
}

INITIAL {
    qt = q10^((celsius-34 (degC))/10 (degC))
	rates(v, cai)
	ab = abinf
}




FUNCTION shiftab(cai (mM))  {
	shiftab = 25 - 55.7 + 136.9*exp(-.28*cai*1e3)
}


FUNCTION peakab(cai (mM))  {
	peakab = 13.7 + 234*exp(-.72*cai*1e3)
}




FUNCTION taufunc(v (mV)) {
	 taufunc = 1 / (          (10*(exp(-v/63.6) + exp (-(150-v)/63.6)))  - 5.2                  )
	 if (taufunc <= 0.2) {	  
	    taufunc = 0.2
	 }

}

PROCEDURE rates(v (mV), cai (mM)) {
	  LOCAL range, vv

	  

	  abinf = -56.449 + 104.52*exp(-.22964*cai*1e3) + 295.68*exp(-2.1571*cai*1e3)

	  abinf = 1/(1+exp((abinf-v)/(25/1.6)))

	  vv = v + 100 - shiftab(cai)
	  abtau = taufunc(vv)
	  range = peakab(cai)-base
	  abtau = ((range*((abtau-.2)/.8)) + base)/qt*2*2

}