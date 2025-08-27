UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (S)  = (siemens)
}

PARAMETER {
	v 		(mV)
	gkbar  = 1.44e-05	(S/cm2) 	

	
        vhalfl = -98.92  (mV)    		
        kl = 10.89       (mV)    		

	
        vhalft=67.0828	 (mV)    		
        at=0.00610779	 (/ms)   		
	bt=0.0817741	 (/ms)	 		

	
        celsius         (degC)  		
        q10 = 1.                              	
}



NEURON {
	SUFFIX kir 			
	USEION k READ ek WRITE ik	
        RANGE  ik, gkbar, vhalfl, kl, vhalft, at, bt, q10 
        GLOBAL linf,taul
}


STATE {
        l
}

ASSIGNED {
        ik                              (mA/cm2)
        gk                              (S/cm2)
        ek                              (mV)
        linf      
        taul
}


INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp	
	gk = gkbar*l			
        ik = gk * ( v - ek )		
}



DERIVATIVE states {     
        rate(v)
        l' =  (linf - l)/taul		
}

PROCEDURE rate(v (mV)) { 
        LOCAL qt
	qt=q10^((celsius-33)/10) 	
        linf = 1/(1 + exp((v-vhalfl)/kl))			
 	taul = 1/(qt *(at*exp(-v/vhalft) + bt*exp(v/vhalft) ))	
}