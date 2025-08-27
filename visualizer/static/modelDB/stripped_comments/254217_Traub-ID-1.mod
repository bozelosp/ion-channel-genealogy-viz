UNITS {
 	(mA) = (milliamp)
 	(mV) = (millivolt)
	(S) = (siemens)		
}

NEURON {
	SUFFIX traub
	NONSPECIFIC_CURRENT i
	RANGE iL,iNa,iK
	RANGE eL, eNa, eK
	RANGE gLbar, gNabar, gKbar
	RANGE v_shft
 }
	
	
PARAMETER {
        gNabar = .03 (S/cm2)	
        gKbar = .015 (S/cm2) 	
        gLbar = 0.00014 (S/cm2) 
        eL = -62.0 (mV) 
        eK = -80 (mV)	
        eNa = 90 (mV)	
        totG = 0
		v_shft = 49.2 
}
 
STATE {
        m h n a b
}
 
ASSIGNED {
        v (mV)
        i (mA/cm2)
        cm (uF)
        iL iNa iK(mA/cm2)
        gNa gK (S/cm2)
	    minf hinf ninf 
		mtau (ms) htau (ms) ntau (ms) 
}


BREAKPOINT {
        SOLVE states METHOD cnexp 
        
        
        gNa = gNabar*h*m*m
		iNa = gNa*(v - eNa)
		
        gK = gKbar*n 
		iK = gK*(v - eK)
        
		iL = gLbar*(v - eL) 
		i = iL + iK + iNa
		
		
		totG = gNa + gK + gLbar
			
}
 

INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
}

? states
DERIVATIVE states {  
	rates(v)
	
	m' = (minf-m)/mtau
	h' = (hinf-h)/htau
	
	n' = (ninf-n)/ntau 
}

? rates
DEFINE Q10 3
PROCEDURE rates(v(mV)) {  
	
	
	LOCAL  alpha, beta, sum, vt, Q
	TABLE 	mtau,ntau,htau,minf,ninf,hinf
	FROM -100 TO 70 WITH 1000
	
	
	
	
	
	
	
		
		
		
		
		
		
		
		
		
		vt = v + v_shft 
		Q = Q10^((35 - 32)/ 10)
		
		if(vt == 13.1){alpha = 0.32*4}
		else{alpha = 0.32*(13.1 - vt)/(exp((13.1 - vt)/4) - 1)}
		if(vt == 40.1){beta = 0.28*5}
		else{beta = 0.28*(vt - 40.1)/(exp((vt - 40.1)/5)-1)}
        sum = alpha + beta
		mtau = 1/sum
		mtau = mtau/Q
        minf = alpha/sum

       
		alpha = 0.128*exp((17 - vt)/18)
		beta = 4/(1 + exp((40 - vt)/5))
        sum = alpha + beta
		htau = 1/sum
		htau = htau/Q
        hinf = alpha/sum

    	
    	if(vt == 35.1){ alpha = 0.016*5 }
		else{alpha =0.016*(35.1 - vt)/(exp((35.1 - vt)/5) - 1)}
		beta = 0.25*exp((20 - vt)/40)
		sum = alpha + beta
        ntau = 1/sum
        ntau = ntau/Q
        ninf = alpha/sum
}