UNITS {
	  (mA) = (milliamp)
	  (mV) = (millivolt)
    (molar) = (1/liter)
    (mM) = (millimolar)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER {		
	  v             (mV)
	  celsius = 34	(degC)
	  gcalbar = 0   (mho/cm2)   
	  ki  = 0.001   (mM)  
	  cai = 5.e-5   (mM)        
	  cao = 2       (mM)        
    tfa = 5                   
    eca = 140     (mV)        
}

NEURON {
	  SUFFIX cal
	  USEION ca READ cai,cao WRITE ica
    RANGE gcalbar, gmax, gcal, minf, taum
}

STATE {	m }                      

ASSIGNED {                       
	  ica   (mA/cm2)
    gcal  (mho/cm2)
    gmax  (mho/cm2) 
    minf
    taum  (ms)
}

INITIAL {                        
    rates(v)
    m = minf
	  gcal = gcalbar*m*h2(cai)
}

BREAKPOINT {
	  SOLVE states METHOD cnexp
	  gcal = gcalbar*m*h2(cai) 
	  ica = gcal*ghk(v,cai,cao)
    if (gcal > gmax) {
        gmax = gcal
    }
}

FUNCTION h2(cai(mM)) {
	  h2 = ki/(ki+cai)
}

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
    LOCAL nu,f
    f = KTF(celsius)/2
    nu = v/f
    ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) { 
    KTF = ((25.(mV)/293.15(degC))*(celsius + 273.15(degC)))
}

FUNCTION efun(z) {
	  if (fabs(z) < 1e-4) {
		    efun = 1 - z/2
	  }else{
		    efun = z/(exp(z) - 1)
	  }
}

FUNCTION alpm(v (mV)) (/ms) {
	  alpm = 0.055(/ms/mV)*(-27.01(mV) - v)/(exp((-27.01(mV)-v)/3.8(mV)) - 1)
}


FUNCTION betm(v (mV)) (/ms) {
    betm =0.94(/ms)*exp((-63.01(mV)-v)/17(mV))
}




DERIVATIVE states {     
    rates(v)
    m' = (minf - m)/taum
}

PROCEDURE rates(v (mV)) { 
    LOCAL a
    TABLE taum, minf FROM -150 TO 150 WITH 300 
    a = alpm(v)
    taum = 1/(tfa*(a+betm(v))) 
    minf = a/(a+betm(v))       
}