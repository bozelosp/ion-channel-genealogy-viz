NEURON {
	SUFFIX HCN
	NONSPECIFIC_CURRENT i
    RANGE gbar, g, taul, linf
    GLOBAL vhalfc, e, vhalfl, kl, vhalft, at, bt, q10, qfact, cAMP, qt
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    (S)  = (siemens)
	(molar) = (1/liter)
	(uM) = (micromolar)		
}

PARAMETER {
	v 		(mV)
    e = -41.9	(mV)    
	gbar=.0001 	(S/cm2) 

				
        vhalfl=-75.3	(mV)	
	kl=7.1			
        









 
    vhalft=30.4	 (mV)    
    at=0.00052	 (/ms)   
	bt=0.2151	 (/ms)	 

				
    celsius         (degC)  
	q10=1.			
	qfact = 1 
	
	Ac = 14 (mV)
	kc 	= -0.35 (/uM)
	chalf 	= -0.8 (uM)
	cAMP = 0 (uM)
	
}

ASSIGNED {
	i (mA/cm2)
        linf      
        taul
        g (S/cm2)
		vhalfc (mV)  
		qt
}

INITIAL {
	qt = q10^((celsius-33)/10)
	rate(v)
	l1=linf
	l2=linf
}

STATE { l1 l2 }

BREAKPOINT {
	SOLVE states METHOD cnexp
	g = gbar*(0.8*l1 + 0.2*l2)  
	i = g*(v-e)

}

DERIVATIVE states {     
        rate(v)
        l1' =  (linf - l1)/taul
		l2' = (linf - l2)/(taul*6.4)  
}

PROCEDURE rate(v (mV)) { 
        

	
	
	if (cAMP <= 0) {
		vhalfc = 0
	}else{
		vhalfc = Ac/(1 + exp((log10(cAMP)-chalf)/kc))  
	}
    linf = 1/(1 + exp((v-vhalfl-vhalfc)/kl))
 	taul = 1/(qt * qfact * (at*exp(-v/vhalft) + bt*exp(v/vhalft) ))
}