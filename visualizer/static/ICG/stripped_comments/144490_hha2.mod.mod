NEURON {
	  SUFFIX hha2
	  USEION na READ ena WRITE ina
	  USEION k READ ek WRITE ik
	  NONSPECIFIC_CURRENT il
	  RANGE gnabar, gkbar, gl, el, gna, gk, gmax
	  RANGE ar2, vhalfs
	  RANGE inf, tau
	  RANGE taus
	  RANGE W
	  GLOBAL taumin
}

UNITS {
	  (mA) = (milliamp)
	  (mV) = (millivolt)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER {                     
    a0r = 0.0003 (/ms)
    b0r = 0.0003 (/ms)
    zetar = 12    
	  zetas = 12   
    gmr = 0.2   
	  ar2 = 1.0               
    
	  taumin = 3   (ms)       
    vvs  = 2     (mV)       
    vhalfr = -60 (mV)       
	  W = 0.016    (/mV)      
    
    
    
    gnabar = 0.0   (mho/cm2)  
	  gkbar = 1.0    (mho/cm2)  
	  gl = 0       (mho/cm2)
	  ena = 60     (mV)       
	  ek = -77     (mV)       
	  el = -70.0   (mV)       
	  celsius = 34 (degC)
	  v            (mV)
    gk           (mho/cm2)
    gna          (mho/cm2)
    gmax         (mho/cm2)
}

STATE {                         
	  m h n s
}

ASSIGNED {			
	  ina    (mA/cm2)
	  ik     (mA/cm2)
	  il     (mA/cm2)
	  inf[4]
	  tau[4] (ms)
}

BREAKPOINT {
	  SOLVE states METHOD cnexp
    gna = gnabar*m*m*h*s
	  ina = gna*(v - ena)                
    gk =  gkbar*n*n
	  ik = gk*(v - ek)                   
	  il = gl*(v - el)                    
    if (gna + gk + gl > gmax) {
        gmax = gna + gk + gl
    }
}

INITIAL {			
	  mhn(v) 
    m = inf[0]
    h = inf[1]
    n = inf[2]
	  s=1
    gna = gnabar*m*m*h*s
	  ina = gna*(v - ena)                
    gk =  gkbar*n*n
	  ik = gk*(v - ek)                   
	  il = gl*(v - el)
    gmax = gk + gna + gl
}

DERIVATIVE states {
	  mhn(v)
	  m' = (inf[0] - m)/tau[0]  
	  h' = (inf[1] - h)/tau[1]  
	  n' = (inf[2] - n)/tau[2]  
	  s' = (inf[3] - s)/tau[3]  
}	

FUNCTION varss(v(mV), i) { 
	  if (i==0) {
        varss = 1 / (1 + exp((v + 44(mV))/(-3(mV))))    
 	  }
	  else if (i==1) {
        varss = 1 / (1 + exp((v + 49(mV))/(3.5(mV))))   
	  }
	  else if (i==2) {	
        varss = 1 / (1 + exp((v + 46.3(mV))/(-3(mV)))) 
        
	  } else {
        
		    varss = alpv(v,vhalfr)
    }
}


FUNCTION alpv(v(mV),vh(mV)) {    
    alpv = (1+ar2*exp((v-vh)/vvs))/(1+exp((v-vh)/vvs))
}

FUNCTION alpr(v(mV)) {       
    alpr = exp((1.e-3)*zetar*(v-vhalfr)*FARADAY/(R*(273.16+celsius))) 
}

FUNCTION betr(v(mV)) {       
    betr = exp((1.e-3)*zetar*gmr*(v-vhalfr)*FARADAY/(R*(273.16+celsius))) 
}

FUNCTION vartau(v(mV), i) (ms) { 
	  LOCAL tmp
	  if (i==0) {
	      vartau = 0.05  
	  }
	  else if (i==1) {
        vartau = 1     
	  }
	  else if (i==2) {
        vartau = 3.5   
    } else {
	      tmp = betr(v)/(a0r+b0r*alpr(v)) 
	      if (tmp<taumin) {tmp=taumin}
	      vartau = tmp      
    }
}	

PROCEDURE mhn(v(mV)) {
    TABLE inf, tau  DEPEND celsius FROM -100 TO 100 WITH 200
	  FROM i=0 TO 3 {
		    tau[i] = vartau(v,i)
		    inf[i] = varss(v,i)
	  }
}