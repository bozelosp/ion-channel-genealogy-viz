NEURON {
	SUFFIX carF
	USEION ca  WRITE ica
        RANGE gcabar, m, h, g, p, eca
	RANGE inf, fac, vha, ka, ta, vhi, ki, ti
	RANGE irtype
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {	
    	v               (mV)
    	celsius = 30	(degC)
	dt              (ms)
    	gcabar = 0.351  (mho/cm2) 
	eca = 10	(mV)      
	vha = -14	(mV)	
	ka = -6.7	(1)	
	ta = 3.6	(ms)	
	vhi = -65	(mV)	
	ki = 11.8	(1)	
	ti = 200	(ms)	
}  

STATE {	m h }               

ASSIGNED {                  
	ica             (mA/cm2)
    	inf[2]
	fac[2]
	tau[2]
	irtype
	g                       
	p
	
}

BREAKPOINT {
	SOLVE states METHOD derivimplicit
	ica = gcabar*m*m*m*h*(v - eca)
	irtype= -ica
	}

INITIAL {
    	m = 0                               
	h = 0.5                             
	states()
	ica = gcabar*m*m*m*h*(v - eca)      
    	irtype=-ica 				
	}

DERIVATIVE states {
	mhn(v*1(/mV))
	m' = (inf[0] - m) / tau[0]
	h' = (inf[1] - h) / tau[1]
}

FUNCTION varss(v, i) {
	if (i==0) {
           varss = 1 / (1 + exp((v-vha)/(ka)))	
	}
	else if (i==1) {    
        varss = 1/ (1 + exp((v-vhi)/(ki)))     
	}
}

FUNCTION vartau(v, i) {
	if (i==0) {
           vartau = ta		
        }
	else if (i==1) {
           vartau = ti		
       }
	
}	

PROCEDURE mhn(v) {LOCAL a, b 

	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)

	}
}