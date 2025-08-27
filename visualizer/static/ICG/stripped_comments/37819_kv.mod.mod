INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  SUFFIX kv
  USEION k READ ek WRITE ik
  RANGE  n, i, gk, gmax
  GLOBAL ninf, ntau
  GLOBAL Ra, Rb
  GLOBAL q10, temp, tadj
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (pS) = (picosiemens)
  (um) = (micron)
} 

PARAMETER {
  gmax = 5   	(pS/um2)	
  v 		(mV)
  
  tha  = 25	(mV)		
  qa   = 9	(mV)		
  
  Ra   = 0.02	(/ms)		
  Rb   = 0.002	(/ms)		

  dt		(ms)
  celsius		(degC)
  temp = 23	(degC)		
  q10  = 2.3			

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



DERIVATIVE  states {   
  rates(v)      
  n' =  (ninf-n)/ntau
}

PROCEDURE rates(v) {  
  

  a = trap0(v,tha,Ra,qa)
  b = trap0(v,tha,-Rb,-qa)
  ntau = 1/tadj/(a+b)
  ninf = a/(a+b)
}

FUNCTION trap0(v,th,a,q) {
  if (fabs(v-th) > 1e-6) {
    trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
  } else {
    trap0 = a * q
  }
}