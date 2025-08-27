INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

DEFINE pspMAX 15				

NEURON {
	POINT_PROCESS glu
	POINTER pre				
	RANGE gmax, i, g, delay
	RANGE delta, flag, erev
	NONSPECIFIC_CURRENT i
	GLOBAL tau,PreThresh
}

ASSIGNED {
	pspSTACK[pspMAX]
	psp0
	psp1
	pspN
	pre
	flag
	delta

}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}

PARAMETER {
	gmax=0.005     (umho)
	erev=0		(mV)
	tau=1.0 			
	delay=1.5	(ms)		
	PreThresh=-20	(mV)		
	v		(mV)
        dt		(ms)

}

STATE {
	X
	Y
}

INITIAL { LOCAL i
	X = 0				
	Y = 0				
	flag = 0			
	delta = 0			
	g = 0
	pspN = 0
	psp0 = 0
	psp1 = 0

  FROM i=0 TO pspMAX {
  pspSTACK[i] = 0 
 }
}

ASSIGNED { i (nA)  g (umho)}

BREAKPOINT {			

  SOLVE dstates METHOD cnexp

  g = gmax * Y / (tau * 0.36787944)   
  i = g*(v - erev)
}

DERIVATIVE dstates {
	CheckThresh()        
	CheckTime()
	X' = delta - X / tau
        Y' = X - Y / tau
}

PROCEDURE CheckThresh() { 

  if (flag) { if (pre < PreThresh) {flag = 0} }
  else {
    if (pre > PreThresh) {
      flag = 1

        AddSpike()

    }
  }
}

PROCEDURE AddSpike() { 

  if (pspN < pspMAX) {
    pspN=pspN+1
    if (psp1 >= pspMAX) { psp1 = 0 }
    pspSTACK[psp1] = t+delay-dt		
    psp1=psp1+1
  }
  else { printf("ERROR
}



PROCEDURE CheckTime() {			
delta = 0
  if ( (pspSTACK[psp0] > 0) && (t > pspSTACK[psp0]) ){		
        delta = 1/dt
        pspN=pspN-1
        pspSTACK[psp0] = 0
        psp0 = psp0+1
        if (psp0 >= pspMAX) { psp0 = 0 }
  }
}