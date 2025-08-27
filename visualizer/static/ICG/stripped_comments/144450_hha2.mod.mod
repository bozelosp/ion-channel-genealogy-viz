NEURON {
	SUFFIX hha2
	USEION na READ ena WRITE ina 
	USEION k READ ek WRITE ik
	NONSPECIFIC_CURRENT il
	RANGE gnabar, gkbar, gl, el,ik,ina
	RANGE ar2, vhalfs
	GLOBAL inf, tau, taumin
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
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
        vvh=-58		(mV) 
        gnabar = 0.0   (mho/cm2)  
	gkbar = 1.0    (mho/cm2)  
	gl = 0       (mho/cm2)  
	el = -70.0   (mV)       

}


ASSIGNED {
	celsius      (degC)
	v            (mV)
	ena          (mV)       
	ek           (mV)       
	ina          (mA/cm2)
	ik           (mA/cm2)
	il           (mA/cm2)
	inf[4]
	tau[4]	     (ms)
}

STATE {
	m 
	h 
	n 
	s
}

INITIAL {                       
	rates(v,ar2)
	m = inf[0]
	h = inf[1]
	n = inf[2]
	s = inf[3]
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ina = gnabar*m*m*h*s*(v - ena) 
	ik = gkbar*n*n*(v - ek)        
	il = gl*(v - el)               
}


DERIVATIVE states {
	rates(v,ar2)
	m' = (inf[0]-m)/tau[0]
	h' = (inf[1]-h)/tau[1]
	n' = (inf[2]-n)/tau[2]
	s' = (inf[3]-s)/tau[3]
}

PROCEDURE rates(v(mV),a2) {
	LOCAL tmp, c
	FROM i=0 TO 2 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
	}
	tau[3] = betr(v)/(a0r*(1+alpr(v))) 
	if (tau[3]<taumin) {tau[3]=taumin} 
	c = alpv(v)
	inf[3] = c+a2*(1-c) 
}


FUNCTION varss(v(mV), i) { 
	if (i==0) {
		varss = 1 / (1 + exp((v + 44)/(-3(mV)))) 
	}
	else if (i==1) {
		varss = 1 / (1 + exp((v + 49)/(3.5(mV))))  
	}
	else if (i==2) {	
		varss = 1 / (1 + exp((v + 46.3)/(-3(mV)))) 

	}
}


FUNCTION alpv(v(mV)) {
         alpv = 1/(1+exp((v-vvh)/vvs))
}

FUNCTION alpr(v(mV)) {       
UNITSOFF
  alpr = exp(1.e-3*zetar*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

FUNCTION betr(v(mV)) {       
UNITSOFF
  betr = exp(1.e-3*zetar*gmr*(v-vhalfr)*9.648e4/(8.315*(273.16+celsius))) 
UNITSON
}

FUNCTION vartau(v(mV), i) (ms){ 
	LOCAL tmp
	if (i==0) {
	   vartau = 0.05(ms)      
	}
	else if (i==1) {
           vartau = 1(ms)       
        }
	else if (i==2) {
            vartau = 3.5(ms)      
       	}
	  
	   
	 
       
}