NEURON {
  SUFFIX km
  USEION k READ ek WRITE ik
  RANGE n, gk, gbar
  RANGE ninf, ntau
  GLOBAL Ra, Rb
  GLOBAL q10, temp, tadj, vmin, vmax
}

UNITS {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (pS) = (picosiemens)
  (um) = (micron)
}

PARAMETER {
  gbar = 10 (pS/um2) 


  tha  = -30  (mV)   
  qa   = 9  (mV)     

  Ra   = 0.001 (/ms) 
  Rb   = 0.001 (/ms) 



  temp = 23  (degC)    
  q10  = 2.3      

  vmin = -120  (mV)
  vmax = 100  (mV)
}

ASSIGNED {
  v     (mV)
  celsius (degC)
  a     (/ms)
  b     (/ms)
  ik    (mA/cm2)
  gk    (pS/um2)
  ek    (mV)
  ninf
  ntau (ms)
  tadj
}

STATE { n }

INITIAL {
  
  
  tadj = q10^((celsius - temp)/(10 (degC)))
  trates(v)
  n = ninf
}

BREAKPOINT {
  SOLVE states METHOD cnexp

  gk = gbar*n
  ik = (1e-4) * gk * (v - ek)
}



DERIVATIVE states { 
  trates(v)         
  n' = (ninf-n)/ntau
  
}

PROCEDURE trates(v (mV)) {  
                      
  TABLE ninf, ntau
  DEPEND celsius, temp, Ra, Rb, tha, qa
  FROM vmin TO vmax WITH 199

  rates(v)

        

        
        
}

UNITSOFF
PROCEDURE rates(v (mV)) {



  a = Ra*qa*efun((tha-v)/qa)
  b = Rb*qa*efun((v-tha)/qa)


  tadj = q10^((celsius - temp)/(10 (degC)))

  ntau = 1/tadj/(a+b)
  ninf = a/(a+b)



}

FUNCTION efun(z) {
  if (fabs(z) < 1e-6) {
    efun = 1 - z/2
  }else{
    efun = z/(exp(z) - 1)
  }
}
UNITSON