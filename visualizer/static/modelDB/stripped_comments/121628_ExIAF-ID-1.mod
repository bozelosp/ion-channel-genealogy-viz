UNITS {
	(mV) = (millivolt)
	(mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX ExIAF 
    NONSPECIFIC_CURRENT i
    GLOBAL spikedur, refact, tauAHP, eAHP, gAHPbar
	RANGE Thr, lastspike
	RANGE gPAS, ePAS, gAHP, AHPon, gON, gOFF, eON, eOFF
    
	RANGE SetCa, AvgCa, Ca, tauDCCa, Gain, PlasThr
	RANGE ScaleFactor, Induction
    GLOBAL SCALE
	GLOBAL tstop
}


PARAMETER {
    v
	i		(mA/cm2)

	gPAS = 0.001			(mho/cm2)
	ePAS = -60				(mV)
	ePASconst = -60
	Thr	   = -40	(mv)
	ThrConst = -40

	spikedur = 1	(ms)
	refact   = 2.0	(ms)

	gONconst  = 1 	(mho/cm2)
	gOFFconst = 1  (mho/cm2)
	eON  = 40		(mV)
	eOFF = -60		(mV)

	tauAHP   = 0.1 	(/ms)		
	gAHPbar = 0.00005 (mho/cm2)	
	eAHP    = -90 	(mv)

	SetCa = 1
	tauDCCa = 10                  
    tstop
    SCALE = 1                   

}

ASSIGNED {
	lastspike
	gAHP 		(mho/cm2)
	AHPon				

	gON		(mho/cm2)
	gOFF		(mho/cm2)

    Ca
	AvgCa
    PlasThr
    Gain
    ScaleFactor
	Induction
    
}


INITIAL {
	ePAS = ePASconst
	gAHP = 0
	AHPon	  = -9e4

	gON = 0
	gOFF = 0

	lastspike = -9e4
	
    Ca=0
    Induction = 0
    

}

BREAKPOINT {
	SOLVE update
	i = gPAS*(v-ePAS) + gAHP*(v-eAHP)+gON*(v-eON)+gOFF*(v-eOFF)
}

PROCEDURE update() { LOCAL q, dv

   gON = 0
   gOFF = 0
   q = (t-lastspike) - spikedur 

   if (q>refact) {				
		if (v>Thr) {			
			gON = gONconst		
			lastspike = t
			Ca = Ca+1
		}
	}
	else if ( q < 0 ) {			
		gON=gONconst			
	} 
	else if (v > 0) {				
		gOFF = gOFFconst
		gAHP = gAHP + gAHPbar
		AHPon = t
	}
	gAHP = gAHP - gAHP*tauAHP*dt



    if ( t>(tstop-(2*dt)+dt/2) ) {
        
        
        if (Induction==0 && SCALE!=0) {
            SOLVE INDUCTION
	    }
	}    



}                              
                       

PROCEDURE INDUCTION() {


    AvgCa = AvgCa+(Ca - AvgCa)/tauDCCa
    PlasThr=(AvgCa^2)/SetCa
    Gain = (Ca-PlasThr)
    ScaleFactor = 1+(0.05*(SetCa-AvgCa)/(SetCa+AvgCa))*SCALE
    
    Induction = 1       
    


}