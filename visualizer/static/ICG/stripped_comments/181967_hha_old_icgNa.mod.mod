NEURON {
	SUFFIX hha_old
	USEION na READ ena WRITE ina 
	
	
	RANGE gnabar, gkbar, gl, el
	RANGE ar2, vhalfs
	RANGE inf, fac, tau
	RANGE taus
	RANGE W
	GLOBAL taumin, ena
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {   
        a0r = 0.0003 (ms)
        b0r = 0.0003 (ms)
        zetar = 12    
	zetas = 12   
        gmr = 0.2   
	ar2 = 1.0               
                                
	taumin = 3   (ms)       
        vvs  = 2     (mV)       
        vhalfr = -60 (mV)       
	W = 0.016    (/mV)      



        gnabar = 1.0   (mho/cm2)  
	gkbar = 1.0    (mho/cm2)  
	gl = 1.0       (mho/cm2)
	ena = 60     (mV)       
	ek = -77     (mV)       
	el = -70.0   (mV)       
	celsius = 34 (degC)
	v            (mV)
        dt           (ms)
}

STATE {                         
	m h n s
}

ASSIGNED {			
	ina (mA/cm2)
	
	
	inf[4]
	fac[4]
	tau[4]
}

BREAKPOINT {
	SOLVE states
	ina = gnabar*m*m*h*s*(v - ena) 
	
	
}

INITIAL {                       
	states()
	s=1
	ina = gnabar*m*m*h*s*(v - ena)
	
	
}

PROCEDURE calcg() {
	mhn(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)   
	h = h + fac[1]*(inf[1] - h)   
	n = n + fac[2]*(inf[2] - n)   
	s = s + fac[3]*(inf[3] - s)   
}	

PROCEDURE states() {	
	calcg()
	VERBATIM
	return 0;
	ENDVERBATIM
}


FUNCTION varss(v, i) { 
	if (i==0) {
		varss = 1 / (1 + exp((v + 40)/(-3))) 
	}
	else if (i==1) {
		varss = 1 / (1 + exp((v + 45)/(3)))  
	}
	else if (i==2) {	
		varss = 1 / (1 + exp((v + 42)/(-2))) 

	} else {
                
		varss = alpv(v,vhalfr)
       }
}


FUNCTION alpv(v(mV),vh) {    
  alpv = (1+ar2*exp((v-vh)/vvs))/(1+exp((v-vh)/vvs))
}

FUNCTION alpr(v(mV)) {       
  alpr = exp(1.e-3*zetar*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION betr(v(mV)) {       
  betr = exp(1.e-3*zetar*gmr*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
}

FUNCTION vartau(v, i) { 
	LOCAL tmp
	if (i==0) {
	   vartau = 0.05      
	}
	else if (i==1) {
           vartau = 0.5       
        }
	else if (i==2) {
            vartau = 2.2      
       	} else {
	     tmp = betr(v)/(a0r+b0r*alpr(v)) 
	     if (tmp<taumin) {tmp=taumin}
	VERBATIM
	ENDVERBATIM
	     vartau = tmp      
       }
}	

PROCEDURE mhn(v) {

	FROM i=0 TO 3 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
}