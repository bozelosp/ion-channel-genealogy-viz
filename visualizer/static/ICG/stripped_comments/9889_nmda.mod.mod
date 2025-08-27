INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS NMDA
  POINTER pre
  RANGE C, R, R0, R1, g, gmax, lastrelease, B, spk
  NONSPECIFIC_CURRENT i
  GLOBAL Cmax, Cdur, Alpha, Beta, Erev, Prethresh, Deadtime, Rinf, Rtau
}

INCLUDE "queue.inc"  

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
  (mM) = (milli/liter)
}

PARAMETER {
  mg    = 1.    (mM)            
  Cmax	= 1	(mM)		
  Cdur	= 5.0	(ms)		
  Alpha	= 0.072	(/ms mM)	
  Beta	= 0.0066 (/ms)		
  Erev	= 0	(mV)		
  Prethresh = 0 		
  Deadtime = 1	(ms)		
  gmax	= 0.001 (umho)		
  vmin = -120     (mV)
  vmax = 100      (mV)
}


ASSIGNED {
  v		(mV)		
  i 		(nA)		
  g 		(umho)		
  C		(mM)		
  B                             
  R				
  R0				
  R1				
  Rinf				
  Rtau		(ms)		
  pre 				
  spk                           
  lastrelease	(ms)		
}

INITIAL {
  initq()                       

  R = 0
  C = 0
  R0 = 0
  R1 = 0
  Rinf = Cmax*Alpha / (Cmax*Alpha + Beta)
  Rtau = 1 / ((Alpha * Cmax) + Beta)
  lastrelease = -9e9
  rates(v)
}

BREAKPOINT {
  rates(v)
  SOLVE release
  g = gmax * R * B
  i = g*(v - Erev)
}

PROCEDURE release() { LOCAL q
  

  if (! spk && pre > Prethresh) { 
    spk = 1
    pushq(t+delay) } 

  if (spk && pre < Prethresh) { 
    spk = 0 }

  q = ((t - lastrelease) - Cdur) 

  
  if (q > Deadtime) {

    if (t >= queu[head]) {      
      popq()                    
      C = Cmax			
      R0 = R
      lastrelease = t
    }
    
  } else if (q < 0) {		

    if (t > queu[head]) { popq() } 

  } else if (C == Cmax) {	
    R1 = R
    C = 0.
  }
  
  if (C > 0) {			
    R = Rinf + (R0 - Rinf) * exptable (- (t - lastrelease) / Rtau)
  } else {			
    R = R1 * exptable (- Beta * (t - (lastrelease + Cdur)))
  }
  
  VERBATIM
  return 0;
  ENDVERBATIM
}

FUNCTION exptable(x) { 
  TABLE  FROM -10 TO 10 WITH 2000
  
  if ((x > -10) && (x < 10)) {
    exptable = exp(x)
  } else {
    exptable = 0.
  }
}

PROCEDURE rates(v(mV)) {
  TABLE B
  DEPEND mg
  FROM vmin TO vmax WITH 200
  
  
  B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}