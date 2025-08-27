NEURON {
	  SUFFIX km
	  USEION k READ ek WRITE ik
	  RANGE n, gk, gbar, gmax
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
	  v 		        (mV)
	  gbar = 10   	(pS/um2)              
	  tha  = -30	  (mV)                  
	  qa   = 9	    (mV)                  
	  Ra   = 0.001	(/ms/mV)              
	  Rb   = 0.001	(/ms/mV)              
	  celsius		    (degC)
	  temp = 23	    (degC)                
	  q10  = 2.3                          
	  vmin = -120	  (mV)
	  vmax = 100	  (mV)
} 


ASSIGNED {
	  a		 (/ms)
	  b		 (/ms)
	  ik 	 (mA/cm2)
	  gk	 (pS/um2)
	  ek   (mV)
	  ninf
	  ntau (ms)	
	  tadj
    gmax (mho/cm2)
}


STATE { n }

INITIAL { 
	  trates(v)
	  n = ninf
    gk = tadj*gbar*n
    gmax = (1e-4) * gk
}

BREAKPOINT {
    SOLVE states METHOD cnexp
	  gk = tadj*gbar*n
	  ik = (1e-4) * gk * (v - ek)
    if (((1e-4) * gk) > gmax) {
        gmax = (1e-4) * gk
    }
} 

LOCAL nexp

DERIVATIVE states {   
    trates(v)     
    n' =  (ninf-n)/ntau*tadj
}

PROCEDURE trates(v (mV)) {  
    
    TABLE ninf, ntau
	  DEPEND celsius, temp, Ra, Rb, tha, qa
	  FROM vmin TO vmax WITH 199
    
	  rates(v)
    tadj = q10^((celsius - temp)/10(degC))  
}


PROCEDURE rates(v (mV)) {  
    
    
    a =  Ra * (v - tha) / (1 - exp(-(v - tha)/qa))
    b = -Rb * (v - tha) / (1 - exp((v - tha)/qa))
    ntau = 1/(a+b)
	  ninf = a*ntau
}