NEURON {
	SUFFIX aabBK
	USEION ca READ ica
	USEION k READ ek WRITE ik
	RANGE gakbar, gabkbar, cal, gak, gabk, atau, ainf, a, abinf, abtau, base, tau
}

UNITS {
	(molar) = (1/liter)
	(mM) = (millimolar)
	(mV) = (millivolt)
	(mA) = (milliamp)
	(S) = (siemens)
	B = .26 (mM-cm2/mA-ms)
}

PARAMETER {
	gakbar = .01	(S/cm2)	
	gabkbar = .01	(S/cm2)	
	ca0 = .00007	(mM)    
	tau = 2		(ms)	
	cal = 0 (mM)
	cascale	= 3		
	base = 1  	(mV)	
}

ASSIGNED {
	v		(mV)
	ek		(mV)
	ik		(mA/cm2)
	ica		(mA/cm2)
	area		(microm2)
      gak		(S/cm2)
      gabk		(S/cm2)
      ainf		(mV)
      atau		(ms)
      abinf		(mV)
      abtau		(ms)
}

STATE { ca_i (mM) a ab }

BREAKPOINT {
	SOLVE state METHOD cnexp
	gak = gakbar*a
	gabk = gabkbar*ab
	ik = (gak+gabk)*(v - ek)
	cal = ca_i
}

DERIVATIVE state {	
	ca_i' = -B*(ica*cascale)-(ca_i-ca0)/tau	      
	rates(v, ca_i)				      
	a' = (ainf-a)/atau	
	ab' = (abinf-ab)/abtau
}

INITIAL {
	ca_i = ca0
	rates(v, ca_i)
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
	  LOCAL range, vv

	  

	  ainf =  -32 + (59.2*exp(-.09*c*1e3)) + (96.7*exp(-.47*c*1e3))
	  ainf = 1/(1+exp((ainf-v)/(25/1.6)))

	  vv = v + 100 - shifta(c)
	  atau = taufunc(vv)
	  range = peaka(c)-1
	  atau = (range*((atau-.2)/.8)) + 1

	  

	  abinf = -56.449 + 104.52*exp(-.22964*c*1e3) + 295.68*exp(-2.1571*c*1e3)

	  abinf = 1/(1+exp((abinf-v)/(25/1.6)))

	  vv = v + 100 - shiftab(c)
	  abtau = taufunc(vv)
	  range = peakab(c)-base
	  abtau = (range*((abtau-.2)/.8)) + base		

}