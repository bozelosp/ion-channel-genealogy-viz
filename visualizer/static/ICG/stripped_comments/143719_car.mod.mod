NEURON {
	SUFFIX car
	USEION ca READ eca WRITE ica

        RANGE gcabar, m, h,ica
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


ASSIGNED {               
        eca (mV)
	ica (mA/cm2)

        inf[2]
	tau[2]		(ms)
        v               (mV)
        celsius 	(degC)
	   
	
	
}


PARAMETER {              
        gcabar = 1.0      (mho/cm2) 
        

       
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

BREAKPOINT {
	SOLVE states METHOD cnexp
	
	ica = gcabar*m*m*m*h*(v - eca)

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
	    varss = 1 / (1 + exp((v+48.5)/(-3(mV)))) 
	}
	else if (i==1) {
             varss = 1/ (1 + exp((v+53)/(1(mV))))    
	}
}

FUNCTION vartau(v(mV), i) (ms){
	if (i==0) {
         vartau = 50(ms)  
   
        }
	else if (i==1) {
          vartau = 5(ms)   
     
       }
	
}