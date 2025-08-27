INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON { 
    NONSPECIFIC_CURRENT i 
	POINT_PROCESS EPlasSyn
    POINTER ampa, nmda
	RANGE B, precell 
	RANGE gAMPA, gmaxAMPA, gNMDA, gmaxNMDA, Erev_1
	GLOBAL Erev_2 
    GLOBAL AMPANMDARATIO, AMPANMDARATIO, terror
    
    
    POINTER ScaleFactor, Induction, lastprespike 
    RANGE scale    
	
    POINTER lastpostspike, Gain
	RANGE LTP, Ca, stdp 	 
    GLOBAL LTPGain, LTDGain, tauLTD, tauLTP

	
} 
 
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}


PARAMETER {

	precell 
 
	Erev_1 = 0		(mV)		
	Erev_2 = 50		(mV)		
    AMPAMAX = 0.1               
    AMPANMDARATIO = 0.1
    terror
    
    
    scale = 1               
    
    stdp = 1
    tauLTD = 10             
    tauLTP = 5              
    LTPGain = 100
    LTDGain = 0.1
    
}


ASSIGNED { 
 
	dt		 (ms)
	v		 (mV)		
	i 		 (nA)		
	C		 (mM)		
	gAMPA	 (umho)	
	gNMDA 	 (umho)	
	gmaxAMPA (umho)	
	gmaxNMDA (umho)	
	B				
		 
	ampa                
    nmda                
    
	
    lastprespike        
	ScaleFactor         
    Induction           

   	
    lastpostspike       
	LTP 
    Ca
	Gain                

}

INITIAL {
    terror = dt/10
	LTP = 0.0
	Ca = 0.0
}

BREAKPOINT {
    SOLVE update
	i = gAMPA*(v-Erev_1) + gNMDA*(v-Erev_2)
}

PROCEDURE update() {

    gAMPA = gmaxAMPA * ampa 
	B = mgblock(v)
	gNMDA = gmaxNMDA * nmda * B 

    if ( Induction==1 && scale==1 && lastprespike>0) {SOLVE SCALE}
    if ( stdp==1 ) {
        if ((lastprespike>t-terror && lastprespike<t+terror) ||  (lastpostspike>t-terror && lastpostspike<t+terror)) {SOLVE STDP}
    }
    VERBATIM
	    return 0;
	ENDVERBATIM
} 


PROCEDURE STDP() { LOCAL dW, dW2
					

	dW = STDPFunc(lastpostspike-lastprespike)
	if (dW>0) {
        dW2=dW*gmaxAMPA*LTPGain*fabs(AMPAMAX-gmaxAMPA)	
	} else {
        dW2=dW*gmaxAMPA*LTDGain		
	}

    
	
	gmaxAMPA = gmaxAMPA + dW2
	
    if (gmaxAMPA>AMPAMAX) {
        gmaxAMPA = AMPAMAX
    }
    gmaxNMDA = gmaxAMPA*AMPANMDARATIO


	VERBATIM
	return 0;
	ENDVERBATIM

}


PROCEDURE SCALE() { LOCAL dummy
    
    gmaxAMPA = gmaxAMPA*ScaleFactor
    if (gmaxAMPA>AMPAMAX) {
        gmaxAMPA = AMPAMAX
    }
    gmaxNMDA = gmaxAMPA*AMPANMDARATIO
    
    

    VERBATIM
    return 0;
    ENDVERBATIM

}


 
FUNCTION mgblock(v(mV)) { 

	TABLE  
	FROM -140 TO 80 WITH 1000 
	if (v>-59) { 
	   mgblock = 1/( 1+exp( (-35-v)/6 ) ) 
	} else { 
	   mgblock = 0 
	} 
} 
 

FUNCTION STDPFunc(x) { 
	TABLE  FROM -100 TO 100 WITH 2000
 
	if (x>-0.9999 && x<0.9999) {
	    STDPFunc=0
	} else if ( x>-100 && x<0)  {
	    STDPFunc = -exp(x/tauLTD)
	} else if (x>0 && x<100) {
	    STDPFunc = exp(-x/tauLTP)
	} else {
		STDPFunc = 0.
	}
}