NEURON {
	  SUFFIX car
	  USEION ca READ eca WRITE ica
    RANGE gcabar, m, h, g, gmax
	  RANGE inf, tau
}

UNITS {
	  (mA) = (milliamp)
	  (mV) = (millivolt)
}

PARAMETER {              
    v             (mV)
    celsius = 34	(degC)
    gcabar = 1.0    (mho/cm2) 
	  eca = 140     (mV)      
}  

STATE {	m h }            

ASSIGNED {               
	  ica    (mA/cm2)
    inf[2]
	  tau[2] (ms)
    g      (mho/cm2)
    gmax   (mho/cm2)
}

BREAKPOINT {
	  SOLVE states METHOD cnexp
    g = gcabar*m*m*m*h
	  ica = g*(v - eca)
    if (g > gmax) {
        gmax = g
    }
}

INITIAL {
    mhn(v)
    m = inf[0]
    h = inf[1]
    g = gcabar*m*m*m*h
    ica = g*(v - eca) 
    gmax = g
}

DERIVATIVE states {
	  mhn(v)
	  m' =  (inf[0] - m)/tau[0]
	  h' =  (inf[1] - h)/tau[1]
}	

FUNCTION varss(v (mV), i) {
	  if (i==0) {
	      varss = 1 / (1 + exp((v+48.5(mV))/(-3(mV)))) 
	  }
	  else if (i==1) {
        varss = 1/ (1 + exp((v+53(mV))/(1(mV))))    
	  }
}

FUNCTION vartau(v (mV), i) (ms) {
	  if (i==0) {
        vartau = 50  
    }
	  else if (i==1) {
        vartau = 5   
    }
	  
}	

PROCEDURE mhn(v (mV)) {LOCAL a, b 
    TABLE inf, tau DEPEND celsius FROM -100 TO 100 WITH 200
  	FROM i=0 TO 1 {
	      tau[i] = vartau(v,i)
		    inf[i] = varss(v,i)
	  }
}