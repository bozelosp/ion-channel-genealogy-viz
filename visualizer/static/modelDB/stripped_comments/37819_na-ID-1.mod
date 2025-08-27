INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  SUFFIX naz
  USEION na READ ena WRITE ina
  RANGE m, h, gna, gmax, i
  GLOBAL tha, thi1, thi2, qa, qi, qinf, thinf
  GLOBAL minf, hinf, mtau, htau
  GLOBAL Ra, Rb, Rd, Rg
  GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
  gmax = 1000   	(pS/um2)	
  vshift = -10	(mV)		
  
  tha  = -35	(mV)		
  qa   = 9	(mV)		
  Ra   = 0.182	(/ms)		
  Rb   = 0.124	(/ms)		

  thi1  = -50	(mV)		
  thi2  = -75	(mV)		
  qi   = 5	(mV)	        
  thinf  = -65	(mV)		
  qinf  = 6.2	(mV)		
  Rg   = 0.0091	(/ms)		
  Rd   = 0.024	(/ms)		

  temp = 23	(degC)		
  q10  = 2.3			

  v 		(mV)
  dt		(ms)
  celsius		(degC)
  vmin = -120	(mV)
  vmax = 100	(mV)
}


UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (pS) = (picosiemens)
  (um) = (micron)
} 

ASSIGNED {
  ina 		(mA/cm2)
  i 		(mA/cm2)
  gna		(pS/um2)
  ena		(mV)
  minf 		hinf
  mtau (ms)	htau (ms)
  tadj
}


STATE { m h }

INITIAL { 
  tadj = q10^((celsius - temp)/10)
  rates(v+vshift)
  m = minf
  h = hinf
}

BREAKPOINT {
  SOLVE states METHOD cnexp
  gna = tadj*gmax*m*m*m*h
  i = (1e-4) * gna * (v - ena)
  ina = i
} 

LOCAL mexp, hexp 

DERIVATIVE states {   
  rates(v+vshift)      
  m' =  (minf-m)/mtau
  h' =  (hinf-h)/htau
}

PROCEDURE rates(vm) {  
  LOCAL  a, b

  a = trap0(vm,tha,Ra,qa)
  b = trap0(-vm,-tha,Rb,qa)

  mtau = 1/tadj/(a+b)
  minf = a/(a+b)

  

  a = trap0(vm,thi1,Rd,qi)
  b = trap0(vm,thi2,-Rg,-qi)
  htau = 1/tadj/(a+b)
  hinf = 1/(1+exp((vm-thinf)/qinf))
}


FUNCTION trap0(v,th,a,q) {
  if (fabs(v-th) > 1e-6) {
    trap0 = a * (v - th) / (1 - exp(-(v - th)/q))
  } else {
    trap0 = a * q
  }
}