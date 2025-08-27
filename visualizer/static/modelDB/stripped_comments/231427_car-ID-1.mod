NEURON {
	SUFFIX car
	USEION ca  WRITE ica
    RANGE gcabar, m, h, g, p, eca
	RANGE inf, fac, tau, k
	GLOBAL irtype
	EXTERNAL Area_canmda
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
	eca = 10		(mV)      

	Area            (cm2)
	k = 1e-06		(mA/nA)

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
	SOLVE states
	ica = gcabar*m*m*m*h*(v - eca)
	irtype= -gcabar*m*m*m*h*(v - eca)
	g = gcabar*m*m*m*h*Area*1e6	
	p = m*m*m*h
	}

INITIAL {
	Area = Area_canmda
    m = 0                               
	h = 0.5                             
	states()
	ica = gcabar*m*m*m*h*(v - eca)      
    irtype=-gcabar*m*m*m*h*(v - eca) 	
	g = gcabar*m*m*m*h*Area*1e6 		
	p = m*m*m*h
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
           varss = 1 / (1 + exp((v+14)/(-6.7)))	
	}
	else if (i==1) {    
        varss = 1/ (1 + exp((v+65)/(11.8)))     
	}
}

FUNCTION vartau(v, i) {
	if (i==0) {
           vartau = 3.6		
        }
	else if (i==1) {
           vartau = 200		
       }
	
}	

PROCEDURE mhn(v) {LOCAL a, b 

	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v,i)
		fac[i] = (1 - exp(-dt/tau[i]))
	}
}