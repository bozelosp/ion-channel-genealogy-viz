INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  SUFFIX kca
  USEION k READ ek WRITE ik
  USEION ca READ cai
  RANGE i, n, gk, gmax
  RANGE ninf, ntau
  GLOBAL Ra, Rb, caix
  GLOBAL q10, temp, tadj, vmin, vmax
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
  cai  		(mM)
  caix = 1	
  
  Ra   = 0.01	(/ms)		
  Rb   = 0.02	(/ms)		

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
  ntau 		(ms)	
  tadj
}


STATE { n }

INITIAL { 
  tadj = q10^((celsius - temp)/10)
  rates(cai)
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
  rates(cai)      
  n' =  (ninf-n)/ntau

}

PROCEDURE rates(cai(mM)) {  
  a = Ra * cai 
  b = Rb
  ntau = 1/tadj/(a+b)
  ninf = a/(a+b)
  
  
}