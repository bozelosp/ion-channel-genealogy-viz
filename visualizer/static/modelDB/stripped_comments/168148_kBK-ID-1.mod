NEURON {
	SUFFIX kBK
	USEION k READ ek WRITE ik
	USEION ca READ cai
	RANGE gpeak, gkact, caPh, caPk, caPmax, caPmin
	RANGE caVhh, CaVhk, caVhmax, caVhmin, k, tau
}


UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) 	= (millimolar)
}


PARAMETER {
		
	gpeak   = 268e-4	(mho/cm2) <0, 1e9>
	
	                                    
	caPh    = 2e-3     (mM)             
	caPk    = 1                         
	caPmax  = 1                         
	caPmin  = 0                         
		
	                                    
	caVhh   = 2e-3    (mM)              
	caVhk   = -0.94208                  
	caVhmax = 155.67 (mV)               
	caVhmin = -46.08 (mV)               
	
	                                    
	                                    
	k       = 17	(mV)
	
	                                    
	                                    
	                                    
	tau     = 1 (ms) <1e-12, 1e9>
	scale   = 100                       
        		

} 	


ASSIGNED {
	v 		(mV)
	ek		(mV)
	ik 		(mA/cm2)
    	cai  		(mM)
	caiScaled	(mM)
	pinf		(1)
}


STATE {
        p
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	ik = gpeak*p* (v - ek)
}

DERIVATIVE states {     
        rate(v, cai)
        p' =  (pinf - p)/tau
}

INITIAL {     
        rate(v, cai)
        p = pinf
}

PROCEDURE rate(v(mV), ca(mM))  {
        caiScaled = ca*scale	
        pinf = P0ca(caiScaled) / ( 1 + exp( (Vhca(caiScaled)-v)/k ) )
}

FUNCTION P0ca(ca(mM)) (1) {
		
	if (ca < 1E-18) { 		
	P0ca = caPmin
	} else {
	P0ca = caPmin + ( (caPmax - caPmin) / ( 1 + (caPh/ca)^caPk ))
	}
}

FUNCTION Vhca(ca(mM)) (mV) {
		
	if (ca < 1E-18) {		
	Vhca = caVhmax
	} else {
	Vhca = caVhmin + ( (caVhmax - caVhmin ) / ( 1 + ((caVhh/ca)^caVhk)) )
	}
}