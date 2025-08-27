UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
   SUFFIX ExIAF
   NONSPECIFIC_CURRENT i
   GLOBAL spikedur, refact, tauAHP, eAHP, gAHPbar, ThrConst
        RANGE Thr, lastspike
        RANGE gPAS, ePAS, gAHP, AHPon, gON, gOFF, eON, eOFF
   
        RANGE SetCa, AvgCa, Ca, tauDCCa, Gain
        RANGE ScaleFactor, Induction
   GLOBAL SCALE, GainConst
        GLOBAL tstop, terror
        RANGE B                          
}


PARAMETER {
   v
   i               (mA/cm2)

   gPAS = 0.0001			
   ePAS = -60              
   ePASconst = -60
   ThrConst = -45			
   Thr  = -45				

   spikedur = 0.6			
   refact   = 3			

   gONconst  = 1		
   gOFFconst = 1		
   eON  = 40			
   eOFF = -53			

   tauAHP   = 0.2			
   gAHPbar = 0.00035		
   eAHP    = -90			

   SetCa = 1
   tauDCCa = 10					
   SCALE = 0						
   GainConst = 0.1					
   tstop
   terror

}

ASSIGNED {
        lastspike
        gAHP            (mho/cm2)
        AHPon                           

        gON             (mho/cm2)
        gOFF            (mho/cm2)

   Ca
   AvgCa
   ScaleFactor
   Induction
   Gain
   B

}


INITIAL {
   ePAS = ePASconst
   gAHP = 0
   AHPon     = -9e4

   gON = 0
   gOFF = 0

   terror = dt/10
   lastspike = -9e4

   Ca=0
   Induction = 0

}

BREAKPOINT {
        SOLVE update
        i = gPAS*(v-ePAS) + gAHP*(v-eAHP)+gON*(v-eON)+gOFF*(v-eOFF)
        B=mgblock(v)
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



    if ( t>(tstop-(2*dt)+terror) ) {
        
        
        if (Induction==0) {
            SOLVE INDUCTION
            }
        }



}                              


PROCEDURE INDUCTION() {


    AvgCa = AvgCa+(Ca - AvgCa)/tauDCCa
    
    
    ScaleFactor = GainConst*(SetCa-AvgCa)*SCALE

    
    Induction = 1
    


}


FUNCTION mgblock(v(mV)) {

        TABLE
        FROM -140 TO 80 WITH 1000
        if (v>-58) {
           mgblock = 1/( 1+exp( (-35-v)/6 ) )
           
        } else {
           mgblock = 0
        }
}