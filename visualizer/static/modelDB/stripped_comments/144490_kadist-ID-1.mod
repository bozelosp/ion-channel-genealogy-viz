UNITS {
	  (mA) = (milliamp)
	  (mV) = (millivolt)
    FARADAY = (faraday) (coulomb)
    R = (k-mole) (joule/degC)
}

PARAMETER { 
	  v               (mV)
    ek = -80        (mV) 
	  celsius = 24	  (degC)
    
	  gkabar = 0      (mho/cm2)  
    vhalfn = -1     (mV)       
    vhalfl = -56    (mV)       
    a0n = 0.1       (/ms)      
    zetan = -1.8    (1)        
    zetal = 3       (1)        
    gmn   = 0.39    (1)        
    gml   = 1       (1)
	  lmin  = 2       (ms)
	  nmin  = 0.1     (ms)
	  pw    = -1      (1)
	  tq    = -40     (mV)
	  qq    = 5       (mV)
	  q10   = 5                  
}

NEURON {
	  SUFFIX kad
	  USEION k READ ek WRITE ik
    RANGE gkabar,gka, gmax
    GLOBAL ninf,linf,taul,taun,lmin
}

STATE {       
	  n l
}

ASSIGNED {       
	  ik   (mA/cm2)
    ninf
    linf      
    taul (ms)
    taun (ms)
    gka  (mho/cm2)
    gmax (mho/cm2)
}

INITIAL {		
	  rates(v)
	  n = ninf
	  l = linf
	  gka = gkabar*n*l
	  ik = gka*(v-ek)	
    gmax = gka
}

BREAKPOINT {
	  SOLVE states  METHOD cnexp
	  gka = gkabar*n*l
	  ik = gka*(v-ek)
    if (gka > gmax) {
        gmax = gka
    }
}

FUNCTION alpn(v(mV)) { LOCAL zeta
    zeta = zetan+pw/(1+exp((v-tq)/qq))
    alpn = exp((1.e-3)*zeta*(v-vhalfn)*FARADAY/(R*(273.16(degC)+celsius))) 
}

FUNCTION betn(v(mV)) { LOCAL zeta
    zeta = zetan+pw/(1+exp((v-tq)/qq))
    betn = exp((1.e-3)*zeta*gmn*(v-vhalfn)*FARADAY/(R*(273.16(degC)+celsius))) 
}

FUNCTION alpl(v(mV)) {
    alpl = exp((1.e-3)*zetal*(v-vhalfl)*FARADAY/(R*(273.16(degC)+celsius))) 
}

FUNCTION betl(v(mV)) {
    betl = exp((1.e-3)*zetal*gml*(v-vhalfl)*FARADAY/(R*(273.16(degC)+celsius))) 
}




DERIVATIVE states {     
    rates(v)
    n' = (ninf - n)/taun
    l' = (linf - l)/taul
}

PROCEDURE rates(v (mV)) {                  
    LOCAL a,qt
    TABLE ninf, taun, linf, taul  DEPEND celsius, vhalfn, vhalfl, a0n, zetan,  zetal, gmn, gml, lmin, nmin,	pw,	tq,	qq,	q10 FROM -100 TO 100 WITH 200
    qt = q10^((celsius-24(degC))/10(degC))         
    a = alpn(v)
    ninf = 1/(1 + a)                   
    taun = betn(v)/(qt*a0n*(1+a))      
	  if (taun<nmin) {taun=nmin}         
    a = alpl(v)
    linf = 1/(1+ a)                    
	  taul = 0.26(ms/mV)*(v+50(mV))                 
	  if (taul<lmin) {taul=lmin}         
}