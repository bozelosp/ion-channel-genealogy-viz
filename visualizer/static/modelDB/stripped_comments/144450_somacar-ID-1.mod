NEURON {
	SUFFIX somacar
	USEION ca READ cai, cao WRITE ica
        RANGE gcabar, m, h, ica

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


PARAMETER {     
        gcabar = 0      (mho/cm2) 
}


ASSIGNED {      
        v               (mV)
 	celsius         (degC)
	ica             (mA/cm2)
	ecar             (mV)      
        inf[2]
	fac[2]
	tau[2]
        cai		(mM)
        cao		(mM)
}

STATE {	
	m 
	h 
} 


INITIAL {
	m = 0    
	h = 1    
	rates(v)

}

BREAKPOINT {

	SOLVE states METHOD cnexp
        ecar = (1e3) * (R*(celsius+273.15))/(2*FARADAY) * log (cao/cai)
	ica = gcabar*m*m*m*h*(v - ecar)
}


DERIVATIVE states {
        rates(v)
	m' = (inf[0]-m)/tau[0]
	h' = (inf[1]-h)/tau[1]
}


PROCEDURE rates(v(mV)) {
	FROM i=0 TO 1 {
		tau[i] = vartau(i)
		inf[i] = varss(v,i)
	}
}



FUNCTION varss(v(mV), i) {
	if (i==0) {
	   varss = 1 / (1 + exp((v+60)/(-3))) 
	}
	else if (i==1) {
           varss = 1/ (1 + exp((v+62)/(1)))   
	}
}

FUNCTION vartau(i) {
	if (i==0) {
           vartau = 100  
        }
	else if (i==1) {
           vartau = 5    
       }
	
}