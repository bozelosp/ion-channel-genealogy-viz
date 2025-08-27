NEURON {
	SUFFIX car
	USEION ca READ eca WRITE ica
        RANGE gcabar, m, h
	RANGE inf, fac, tau
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

PARAMETER {              
        v               (mV)
        celsius = 34	(degC)
	dt              (ms)
        gcabar = 0      (mho/cm2) 
	eca = 140       (mV)      
        }  

STATE {	m h }            

ASSIGNED {               
	ica (mA/cm2)
        inf[2]
	fac[2]
	tau[2]
}

BREAKPOINT {
	SOLVE states
	ica = gcabar*m*m*m*h*(v - eca)
	}

INITIAL {
        m = 0    
	h = 1    
	states()
	ica = gcabar*m*m*m*h*(v - eca) 
        }

PROCEDURE calcg() {
	mhn(v*1(/mV))
	m = m + fac[0]*(inf[0] - m)
	h = h + fac[1]*(inf[1] - h)
	}	

PROCEDURE states() {	
	calcg()
	VERBATIM
	return 0;
	ENDVERBATIM
}

FUNCTION varss(v, i) {
	if (i==0) {
	    varss = 1 / (1 + exp((v+48.5)/(-3))) 
	}
	else if (i==1) {
             varss = 1/ (1 + exp((v+53)/(1)))    
	}
}

FUNCTION vartau(v, i) {
	if (i==0) {
           vartau = 50  
        }
	else if (i==1) {
           vartau = 5   
       }
	
}	

PROCEDURE mhn(v) {LOCAL a, b 

	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
}