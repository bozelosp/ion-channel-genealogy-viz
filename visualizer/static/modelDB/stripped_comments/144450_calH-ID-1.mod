NEURON {
	SUFFIX calH
	USEION ca READ cai, cao WRITE ica
        RANGE gcalbar, m, h, ica
	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) =	(millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


PARAMETER { 
        ki     = 0.1  (mM)            
        gcalbar = 0     (mho/cm2) 
}


ASSIGNED {                        
	ica 		(mA/cm2)
        inf[2]				
	tau[2]    	(ms)		
        v               (mV)
        celsius	        (degC)
	ecan             (mV)            
	cai             (mM)      
	cao             (mM)      
}

STATE {	
	m 
	h 
}                  

INITIAL {
        rates(v)
        m = 0    
	h = 1    
}

FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}



BREAKPOINT {
	SOLVE states METHOD cnexp
	ecan = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcalbar*m*m*m*h*h2(cai)*(v - ecan)       
}

DERIVATIVE states {
	rates(v)
	m' = (inf[0]-m)/tau[0]
	h' = (inf[1]-h)/tau[1]
}


PROCEDURE rates(v(mV)) {LOCAL a, b 
	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
	}
}


FUNCTION varss(v(mV), i) {
	if (i==0) { 

             varss = 1 / (1 + exp((v+37.7)/(-1(mV))))  
	}
	else if (i==1) { 
             varss = 1 / (1 + exp((v+41)/(0.5(mV)))) 
	}
}

FUNCTION vartau(v(mV), i) (ms){
	if (i==0) {

           vartau = 3.5(ms)  
        }
	else if (i==1) {

           vartau = 20(ms)   
        }
}