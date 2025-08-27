UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX KDR
	USEION k READ ki,ek WRITE ik
	RANGE gk
	GLOBAL rest,activate_Q10,Q10,gmaxQ10,rate_k,gmax_k,temp1,temp2,tempb
}

PARAMETER {
        v (mV)
	dt (ms)
	gk    = 3.842910637e-03 (mho/cm2)
	rest  = -60.0 (mV) 
	ek
	ki
	celsius
	
	activate_Q10 = 1
	Q10 = 1.200000603e+00
	gmaxQ10 = 1.200000603e+00
	temp1 = 19.0 (degC)
	temp2 = 29.0 (degC)
	tempb = 23.0 (degC)
}

STATE {
        n 
}

ASSIGNED { 
	ik (mA/cm2)
	alphan (/ms)
	betan (/ms)
	rate_k
	gmax_k
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik   = (gk*gmax_k)*n*(v-ek)
}

UNITSOFF

INITIAL {
	LOCAL ktemp,ktempb,ktemp1,ktemp2
	if (activate_Q10>0) {
          rate_k = Q10^((celsius-tempb)/10)
          gmax_k = gmaxQ10^((celsius-tempb)/10)
	}else{
	  
          
          
	  rate_k = 1.60
	  gmax_k = 1.60
	}
        settables(v)
	n = alphan/(alphan+betan)
}

DERIVATIVE states {
	settables(v)      
	n' = alphan * (1-n) - betan * n
}

PROCEDURE settables(v) {  
                          
                          
        LOCAL vadj
        TABLE alphan, betan DEPEND rest,celsius FROM -100 TO 100 WITH 400
	vadj  = v - rest + 0.60650122

	        
	alphan = rate_k * 0.01 * vtrap((35.1-vadj),5.0)
        betan = rate_k * 0.156 * exp((20.0-vadj)/40.0)
}

FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrap = y*(1 - x/y/2)
        }else{
                vtrap = x/(exp(x/y) - 1)
        }
}
 
UNITSON