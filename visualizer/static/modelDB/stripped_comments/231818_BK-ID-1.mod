NEURON {
	SUFFIX BK
	USEION ca READ cai
	USEION k READ ek WRITE ik
	USEION nca READ ncai VALENCE 0
	RANGE gakbar, gabkbar, gak, gabk, atau, ainf, a, ab, abinf, abtau,  ik, diff, acai
	GLOBAL base
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
}

PARAMETER {
	diff = 1 (1)
	gakbar = .01	(S/cm2)	
	gabkbar = .01	(S/cm2)	
	base = 4  	(mV)	
	
}

ASSIGNED {
	v				(mV)
	ek			(mV)
	
	ik				(mA/cm2)
      gak			(S/cm2)
      gabk		(S/cm2)
      ainf	
      atau		(ms)
      abinf	
      abtau		(ms)
	  cai 			(mM)
	  ncai		(mM)
	  acai		(mM)

}

STATE {ab a} 

BREAKPOINT {
	SOLVE state METHOD cnexp
	gak = gakbar*a
	gabk = gabkbar*ab
	ik = (gabk+gak)*(v - ek)                               
}

DERIVATIVE state {
	
	
	acai = (ncai) / diff 

	if (acai < cai)
		{acai = cai}
		
	rates(v,acai)				      
	a' = (ainf-a)/atau	
	ab' = (abinf-ab)/abtau
}

INITIAL {
	rates(v, cai)
	a = ainf
	ab = abinf
}


FUNCTION shifta(ca (mM))  {
	shifta = 25 - 50.3 + (107.5*exp(-.12*ca*1e3))
}


FUNCTION peaka(ca (mM))  {
	peaka = 2.9 + (6.3*exp(-.36*ca*1e3))
}




FUNCTION shiftab(ca (mM))  {
	shiftab = 25 - 55.7 + 136.9*exp(-.28*ca*1e3)
}


FUNCTION peakab(ca (mM))  {
	peakab = 13.7 + 234*exp(-.72*ca*1e3)
}




FUNCTION taufunc(v (mV)) {
	 taufunc = 1 / (          (10*(exp(-v/63.6) + exp (-(150-v)/63.6)))  - 5.2                  )
	 if (taufunc <= 0.2) {	  
	    taufunc = 0.2
	 }
}

PROCEDURE rates(v (mV), c (mM)) { 
	  LOCAL range, vv,ashift, bshift

	  

	  ashift =  -32 + (59.2*exp(-.09*c*1e3)) + (96.7*exp(-.47*c*1e3))
	  ainf = 1/(1+exp((ashift-v)/(25/1.6)))

	  vv = v + 100 - shifta(c)
	  atau = taufunc(vv)
	  range = peaka(c)-1
	  atau = (range*((atau-.2)/.8)) + 1

	  

	  bshift = -56.449 + 104.52*exp(-.22964*c*1e3) + 295.68*exp(-2.1571*c*1e3)

	  abinf = 1/(1+exp((bshift-v)/(25/1.6)))

	  vv = v + 100 - shiftab(c)
	  abtau = taufunc(vv)
	  range = peakab(c)-base
	  abtau = (range*((abtau-.2)/.8)) + base		

}