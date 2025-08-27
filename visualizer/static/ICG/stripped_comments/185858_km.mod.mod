NEURON {
  SUFFIX km
  USEION k READ ek WRITE ik
  RANGE n, gk, gmax, i
  RANGE ninf, ntau, tadj
  GLOBAL Ra, Rb, ek
  GLOBAL q10, temp, vmin, vmax
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (pS) = (picosiemens)
  (um) = (micron)
} 

PARAMETER {
  gmax = 10   	(pS/um2)	
  v 		(mV)
  
  tha  = -30	(mV)		
  qa   = 9	(mV)		
  
  Ra   = 0.001	(/ms)		
  Rb   = 0.001	(/ms)		

  dt		(ms)
  celsius		(degC)
  temp = 23	(degC)		
  q10  = 2.3			

  vmin = -120	(mV)
  vmax = 100	(mV)
} 


ASSIGNED {
  a		(/ms)
  b		(/ms)
  i 		(mA/cm2)
  ik 		(mA/cm2)
  gk		(pS/um2)
  ek		(mV)
  ninf
  ntau (ms)	
  tadj
}

STATE { n }

INITIAL { 
  tadj = q10^((celsius - temp)/10)
  rates(v)
  n = ninf
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  gk = tadj*gmax*n
  i = (1e-4) * gk * (v - ek)
  ik = i
} 

LOCAL nexp

DERIVATIVE states {   
  rates(v)      
  n' = (ninf-n)/ntau

}

PROCEDURE rates(v) {  
  

  a = Ra * (v - tha) / (1 - exp(-(v - tha)/qa))
  b = -Rb * (v - tha) / (1 - exp((v - tha)/qa))

  ntau = 1/tadj/(a+b)
  ninf = a/(a+b)
}