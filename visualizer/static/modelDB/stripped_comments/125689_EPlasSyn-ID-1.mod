INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
   NONSPECIFIC_CURRENT i
        POINT_PROCESS EPlasSyn
   POINTER ampa, nmda
   POINTER PreAvgCa
   POINTER postB
        RANGE precell
        RANGE gAMPA, gmaxAMPA, gNMDA, gmaxNMDA
        GLOBAL Erev_1, Erev_2
   GLOBAL AMPANMDARATIO, AMPANMDARATIO

   
   POINTER ScaleFactor, Induction, lastprespike, AMPAMAX
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
   Erev_1 = 0           (mV)            
        Erev_2 = 50             (mV)            
   AMPANMDARATIO = 0.1

   
   scale = 0               

   
   tauLTP= 20              
   tauLTD= 20              
   gainLTP = 0.05           
   gainLTD = 0.1
   stdp = 0                
   plast = 0                
   terror
   lastanyspike

}


ASSIGNED {

        dt               (ms)
        v                (mV)              
        i                (nA)              
        C                (mM)              
        gAMPA    (umho)    
        gNMDA    (umho) 
        gmaxAMPA (umho) 
        gmaxNMDA (umho) 

   ampa                
   nmda                
   PreAvgCa            
   postB               

   
   lastprespike         
   ScaleFactor          
   Induction            
   AMPAMAX

        
        lastpostspike        


}

INITIAL {
   terror=dt/10
   plast = 0
   lastanyspike=-8e4
}

BREAKPOINT {
   SOLVE update
   i = gAMPA*(v-Erev_1) + gNMDA*(v-Erev_2)
}

PROCEDURE update() {

   gAMPA = gmaxAMPA * ampa
   gNMDA = gmaxNMDA * nmda * postB

   if (stdp==1) {SOLVE STDP}
   if (Induction==1) {SOLVE SCALE}

   VERBATIM
           return 0;
   ENDVERBATIM

}



PROCEDURE SCALE() { LOCAL dummy

   
   
   if (scale>0) {
      gmaxAMPA = gmaxAMPA*( (ScaleFactor)*scale*PreAvgCa+1 )
   }

   
   gmaxAMPA=gmaxAMPA+gmaxAMPA*plast

   
   if (scale > 0) {
		if (gmaxAMPA>AMPAMAX) {
			gmaxAMPA = AMPAMAX
		}
		if (gmaxAMPA<0) {
			gmaxAMPA=0
		}
	}

   
   gmaxNMDA = gmaxAMPA*AMPANMDARATIO

   

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
        DEPEND gainLTP, gainLTD
        FROM -100 TO 100 WITH 2010
        ISI=ISI-3
        if (ISI<=0.1+terror && ISI>-100) {
           STDPFunc = -exp(ISI/tauLTD)*gainLTD
        } else if (ISI>0.1+terror && ISI<100) {
           STDPFunc = exp(-ISI/tauLTP)*gainLTP
        } else {
           STDPFunc = 0
        }
}