INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON { 
   NONSPECIFIC_CURRENT i 
	POINT_PROCESS IPlasSyn
   POINTER gaba
	RANGE precell 
	RANGE gGABA, gmaxGABA
	GLOBAL Erev
	RANGE precell 
   POINTER PreAvgCa
   
   
   POINTER ScaleFactor, Induction, lastprespike, GABAMAX 
   RANGE scale
   
   
   POINTER lastpostspike
   GLOBAL tauLTP, tauLTD, gainLTP, gainLTD
   RANGE stdp, plast, lastanyspike
   GLOBAL terror
} 
 
UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
}


PARAMETER {

	precell 
	Erev = -70		(mV)		
   
   
   scale = 0               

   
   tauLTP= 15              
   tauLTD= 5              
   gainLTP =0.8           
   gainLTD = 1
   stdp = 0                
   plast = 0                
   terror
   lastanyspike

}


ASSIGNED { 
 
	dt		 (ms)
	v		 (mV)		   
	i 		 (nA)		   
	C		 (mM)		   

   gGABA	 (umho)	   
	gmaxGABA (umho)	
	gaba              

   PreAvgCa            
   
	
   lastprespike         
	ScaleFactor          
   Induction            
   GABAMAX

 	
 	lastpostspike        

   
}

INITIAL {
   terror=dt/10 
   plast = 0
   lastanyspike=-8e4
}

BREAKPOINT {
   SOLVE update
	i = gGABA*(v-Erev)
}

PROCEDURE update() {

	gGABA = gmaxGABA * gaba 
 
   if (stdp==1) {SOLVE STDP}
	if (Induction==1) {SOLVE SCALE}
	
   VERBATIM
	   return 0;
	ENDVERBATIM

} 



PROCEDURE SCALE() { LOCAL dummy

   
   
   if (scale>0) {
      
      
      gmaxGABA = gmaxGABA*( (ScaleFactor*(-0.25))*scale*PreAvgCa+1 )
   }
   
   
   gmaxGABA=gmaxGABA+gmaxGABA*plast
   
   
   if (scale > 0) {
		if (gmaxGABA>GABAMAX) {
			gmaxGABA = GABAMAX
		}
   if (gmaxGABA<GABAMAX/1000) {
			gmaxGABA=GABAMAX/1000
		}
	}

   

   VERBATIM
   return 0;
   ENDVERBATIM

}

PROCEDURE STDP() { LOCAL dummy

   if (lastprespike-terror>lastanyspike || lastpostspike-terror>lastanyspike) {
      dummy=lastpostspike-lastprespike
      plast = plast + STDPFunc(dummy)
      
      if (dummy+terror>=0) {
         lastanyspike=lastpostspike
      } else if (dummy<0) {
         lastanyspike=lastprespike
      }
      
   }   
   
   
   VERBATIM
      return 0;
   ENDVERBATIM

}

FUNCTION STDPFunc(ISI) { 
	TABLE  
	FROM -100 TO 100 WITH 2010
	
	ISI=ISI-3
	if (ISI<=0.0+terror && ISI>-100) {
      STDPFunc = exp(ISI/tauLTP)*gainLTP
   } else if (ISI>0.0+terror && ISI<100) {
      STDPFunc = -exp(-ISI/tauLTD)*gainLTD
   } else {
      STDPFunc = 0
   }
	
	
	
   
   
   
   
   
   

   
	
   
   
   
   
   
   

}