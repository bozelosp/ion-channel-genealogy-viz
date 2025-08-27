INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  SUFFIX Nca
  USEION ca READ eca WRITE ica
  RANGE i, m, h, gca, gmax
  RANGE minf, hinf, mtau, htau
  GLOBAL q10, temp, tadj, vmin, vmax, vshift
}

PARAMETER {
  gmax = 0.1   	(pS/um2)	
  vshift = 0	(mV)		

  cao  = 2.5	(mM)	        
  cai		(mM)
  
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
  FARADAY = (faraday) (coulomb)
  R = (k-mole) (joule/degC)
  PI	= (pi) (1)
} 

ASSIGNED {
  i 		(mA/cm2)
  ica 		(mA/cm2)
  gca		(pS/um2)
  eca		(mV)
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
  gca = tadj*gmax*m*m*h
  i = (1e-4) * gca * (v - eca)
  ica = i
} 

LOCAL mexp, hexp


  
  
  
  
  

  DERIVATIVE states {
    rates(v+vshift)      
    m' =  (minf-m)/mtau
    h' =  (hinf-h)/htau
  }

  PROCEDURE rates(vm) {  
    LOCAL  a, b

    a = 0.055*(-27 - vm)/(exp((-27-vm)/3.8) - 1)
    b = 0.94*exp((-75-vm)/17)
    
    mtau = 1/tadj/(a+b)
    minf = a/(a+b)

    

    a = 0.000457*exp((-13-vm)/50)
    b = 0.0065/(exp((-vm-15)/28) + 1)

    htau = 1/tadj/(a+b)
    hinf = a/(a+b)
  }

  FUNCTION efun(z) {
    if (fabs(z) < 1e-4) {
      efun = 1 - z/2
    }else{
      efun = z/(exp(z) - 1)
    }
  }