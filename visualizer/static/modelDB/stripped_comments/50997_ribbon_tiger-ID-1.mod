DEFINE NumSites 10 

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
     POINT_PROCESS ribbon_tiger
     EXTERNAL mod_modulator


     POINTER preCA1, preCA2

     RANGE preCAnow
     RANGE R, RM, RA, RA2, RdA, RdA2, RMA, RMA2, RdMA, RdMA2, O, OM
     RANGE g, gmax, erev, rb, rBb, rDb, rDBb, rMb, rMBb, rm, rBm, rB2m, rDBm, rDB2m
     RANGE Rb, RBu, RBd, RBDr, RBb, RB2u
     RANGE RB2d, RB2Dr, RB2o, RB2c, RDBb, RDB2u
     RANGE Rm, RMum, RBm, RBMum, RB2m, RB2Mum
     RANGE RMb, RMBu, RMBb, RMB2u, RMB2o, RMB2c
     RANGE RMB2d, RMB2Dr, RMBd, RMBDr, RMDBb, RMDB2u
     RANGE RDBMum, RDBm, RDB2Mum, RDB2m
     RANGE AbsRefract,  gluConc_dAMPA, pr, StanddAMPA, rate_constant

     RANGE g2, gmaxN, erevN, rb2, mg
     RANGE C0, C1, C2, D, OO, B
     RANGE Rb2, Ru, Rd, Rr, Ro, Rc
     RANGE gluConc_NMDA, StandNMDA

     RANGE XC,w,A
     RANGE Max_RVP,  rate_constantIN
     RANGE Vmax, k, n
     RANGE UFP, Orate
     RANGE S,W,K1,powr
     GLOBAL total_count



     RANGE xC0, xC1, xC2, xD, xOO

     GLOBAL vmin, vmax
     NONSPECIFIC_CURRENT i
}

UNITS {
     (nA) = (nanoamp)
     (mV) = (millivolt)
     (pS) = (picosiemens)
     (umho) = (micromho)
     (mM) = (milli/liter)
     (uM) = (micro/liter)
}

PARAMETER {
     preCA1 (mM) 
     preCA2 (mM) 
     preCAnow
     erev = 0 (mV)	
     erevN = 0 (mV)     
     gmax = 0 (umho)	
     gmaxN = 0 (umho)   
     StanddAMPA = 0 (mM) 
     StandNMDA = 0.001 (mM)  
	xC0 = 0.442216   
	xC1 = 0.23271    
	xC2 = 0.12248    
	xD =  0.15075    
	xOO = 0.0519     
     mg = 0 (mM)        
        AbsRefract = 0 (ms) 
        rate_constant (1/s) 
        rate_constantIN = 0.25 (1/s) 
        XC 
        w  
        A  
        Max_RVP = 5 
        Vmax = 1840.22 
        k = 86.48823 
        n = 3.30393  
        S
        W
        K1
        powr



     Rb  = 2e-2     (/uM /ms) 
     RBu  = 3e-1    (/ms)     
     RBd  = 1e0     (/ms)     
     RBDr  = 3e-1   (/ms)     
     RBb  =1e-2     (/uM /ms) 
     RB2u  = 1e2    (/ms)     
     RB2d = 8e0     (/ms)     
     RB2Dr = 2e-4   (/ms)     
     RB2o = 3e1     (/ms)     
     RB2c = 1.5e0   (/ms)     
     RDBb = 1e-2    (/uM /ms) 
     RDB2u = 8.3e-3      (/ms)     


     Rm  = 0 (/uM /ms)       
     RMum  = 0        (/ms)   
     RBm = 0         (/uM /ms)       
     RBMum =0         (/ms)   
     RB2m = 0        (/uM /ms)       
     RB2Mum = 0       (/ms)   
     RMb = 0  (/uM /ms)      
     RMBu = 0        (/ms)   
     RMBb = 0        (/uM /ms)       
     RMB2u = 0        (/ms)   
     RMB2o = 0        (/ms)   
     RMB2c = 0     (/ms)   
     RMB2d = 0        (/ms)   
     RMB2Dr = 0      (/ms)   
     RMBd = 0 (/ms)   
     RMBDr = 0       (/ms)   
     RMDBb = 0       (/uM /ms)       
     RMDB2u = 0 (/ms)      
     RDBm = 0        (/uM /ms)       
     RDBMum = 0      (/ms)   
     RDB2m = 0       (/uM /ms)       
     RDB2Mum = 0 (/ms)       


        Rb2     = 5e-3    (/uM /ms)     
        Ru      = 9.5e-3   (/ms)        
        Rd      = 16e-3    (/ms)        
        Rr      = 13e-3    (/ms)        
        Ro      = 25e-3     (/ms)       
        Rc      = 59e-3     (/ms)       

        vmin = -200     (mV)
        vmax = 100      (mV)

}




ASSIGNED {
     gluConc_dAMPA    (mM) 
     gluConc_NMDA     (mM) 
     v		(mV)	
     i		(nA)	
     g		(umho)	
     g2         (umho)  
     mod_modulator (mM)
     rb		(/ms)	
     rb2        (/ms)
     rBb		(/ms)	
     rDBb	(/ms)	
     rMb		(/ms)	
     rMBb	(/ms)	
     rm		(/ms)	
     rBm		(/ms)	
     rB2m	(/ms)	
     rDBm	(/ms)	
     rDB2m	(/ms)	
     dt         (ms)
     last_quanta[NumSites] (ms)   
     release_start[NumSites] (ms) 
     RVP_Size[NumSites] 
     RVP_out[NumSites] 
     total_count 
     Orate      (/s)  
     rate[NumSites] (/s) 
     rateIN[NumSites] (/s) 
     UFP 
}

STATE {

     R		
     RM		
     RA		
     RA2		
     RdA		
     RdA2		
     RMA		
     RMA2		
     RdMA		
     RdMA2		
     O		
     OM		


        C0              
        C1              
        C2              
        D               
        OO              

        B               
}


INITIAL {
     R = 1
FROM i=0 TO NumSites-1{     
           last_quanta[i] = 40
           release_start[i] = 10000 
RVP_Size[i] = Max_RVP
     }
        rates(v)

        C0 = xC0 
        C1 = xC1
        C2 = xC2
        D = xD
        OO = xOO

total_count = 0
UFP = 0
}

BREAKPOINT {
rates(v)

gluConc_dAMPA = StanddAMPA
gluConc_NMDA = StandNMDA

Orate = 0

UFP = 0
FROM i=0 TO NumSites-1
{ 
if (RVP_Size[i] == Max_RVP)
   {UFP = UFP + 1}
}



FROM i=0 TO NumSites-1 
  {

if (RVP_Size[i] < Max_RVP) {preCAnow = preCA2}
else
{preCAnow = preCA1}



if (RVP_Size[i] < 1)

       { S = 2.84985  W = 0.0498303 K1 = 18.3892 powr=4
         rate[i] = S * (W * preCAnow * 1000) / ( ((preCAnow * 1000) / K1) + 1)^powr }

else 
{




  rate_constant =  (Vmax * (1000 * preCAnow)^n) / (k^n + (1000 * preCAnow)^n)

  if (RVP_Size[i] == Max_RVP)
   { rate[i] = UFP * rate_constant }
  else
   { rate[i] =  RVP_Size[i] * rate_constant }
}
 
Orate = Orate + rate[i]


  last_quanta[i] = last_quanta[i] + dt
  if ( ( scop_random() < (rate[i] / (1000/dt) ) ) && (last_quanta[i] >= AbsRefract) )
     {
       release_start[i] = t
       last_quanta[i] = 0
       RVP_out[i] = 1
       total_count = total_count + 1
     }



XC=0.09 w=1.28093 A=0.09262
gluConc_dAMPA= gluConc_dAMPA + transmitter(t - release_start[i])
XC=0.19 w=1.29098 A=0.0256
gluConc_NMDA = gluConc_NMDA + transmitter(t - release_start[i])
  }



FROM i=0 TO NumSites-1
{
RVP_Size[i] = RVP_Size[i]  - RVP_out[i]
if (RVP_Size[i] < 0) {RVP_Size[i] = 0}
if (RVP_Size[i] > Max_RVP){RVP_Size[i] = Max_RVP}

if (RVP_Size[i] < Max_RVP){
rateIN[i] = rate_constantIN  *  (Max_RVP - RVP_Size[i]) 
 if (scop_random() < rateIN[i] / (1000 / dt))
  { RVP_Size[i] = RVP_Size[i] + 1 }
}
RVP_out[i] = 0
}

SOLVE kstates METHOD sparse

     g = gmax *(O+OM) 
     g2 = gmaxN * OO * B 
     i = (g * (v - erev) ) + (g2 * (v - erevN) )
}


KINETIC kstates {
     rb = Rb * (1e3) * gluConc_dAMPA		
     rBb = RBb * (1e3) * gluConc_dAMPA		
     rDBb = RDBb * (1e3) * gluConc_dAMPA		
     rMb = RMb * (1e3) * gluConc_dAMPA 		
     rMBb = RMBb * (1e3) * gluConc_dAMPA 
     rm  = Rm * (1e3) * mod_modulator		
     rBm = RBm * (1e3) * mod_modulator   	
     rB2m = RB2m * (1e3) * mod_modulator	
     rDBm = RDBm * (1e3) * mod_modulator	
     rDB2m = RDB2m * (1e3) * mod_modulator	

     ~ R <-> RA     (rb,RBu)
     ~ RA <-> RdA   (RBd,RBDr)
     ~ RA <-> RA2   (rBb,RB2u)
     ~ RA2 <-> RdA2 (RB2d,RB2Dr)
     ~ RA2 <-> O    (RB2o,RB2c)
     ~ RdA <-> RdA2 (rDBb,RDB2u)
     ~ R <-> RM     (rm,RMum)
     ~ RA <-> RMA   (rBm,RBMum)
     ~ RA2 <-> RMA2 (rB2m,RB2Mum)
     ~ RM <-> RMA   (rMb,RMBu)
     ~ RMA <-> RMA2 (rMBb,RMB2u)
     ~ RMA2 <-> OM   (RMB2o,RMB2c)
     ~ RMA2 <-> RdMA2 (RMB2d,RMB2Dr)
     ~ RMA <-> RdMA (RMBd,RMBDr)
     ~ RdMA <-> RdMA2 (RMDBb,RMDB2u)
     ~ RdA <-> RdMA (rDBm,RDBMum)
     ~ RdA2 <-> RdMA2 (rDB2m,RDB2Mum)

     CONSERVE R+RM+RA+RA2+RdA+RdA2+RMA+RMA2+RdMA+RdMA2+O+OM = 1


       rb2 = Rb2 * (1e3) * gluConc_NMDA

        ~ C0 <-> C1     (rb2,Ru)
        ~ C1 <-> C2     (rb2,Ru)
        ~ C2 <-> D      (Rd,Rr)
        ~ C2 <-> OO     (Ro,Rc)

        CONSERVE C0+C1+C2+D+OO = 1
}

PROCEDURE rates(v(mV)) {
        TABLE B
        DEPEND mg
        FROM vmin TO vmax WITH 200

        

        B = 1 / (1 + exp(0.062 (/mV) * -v) * (mg / 3.57 (mM)))
}


FUNCTION transmitter(x) {

        if (x <= 0 || x > 10) {
             transmitter = 0
        }else{
             transmitter = A/(sqrt(2*3.1415927)*w*x) * exp(-log(x/XC)^2 / (2*w^2))
        }
}