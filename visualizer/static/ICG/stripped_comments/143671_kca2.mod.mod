NEURON {
 	SUFFIX kca2
 	USEION k READ ek WRITE ik
 	USEION ca READ cai 
  
  RANGE n, g,ik,cai,ica,icaL,depth1,taur1,depth2,taur2
 	GLOBAL Ra, Rb, caix
 }
 
 UNITS {
 	(mA) = (milliamp)
 	(mV) = (millivolt)
 	(S) = (siemens)
 	(um) = (micron)
 	(molar) = (1/liter)			
 	(mM)	= (millimolar)
 	(msM)	= (ms mM)
 	FARADAY = (faraday) (coulomb)
 } 
 
 PARAMETER {
 	g = 0.03   	(S/cm2)	
 	v 		(mV)
 	cai  		(mM)
 	caix = 2	
  cainf=0.0001
 	depth1	= .1	(um)		
 	taur1	= 20	(ms)		
 	depth2	= 10	(um)		
 	taur2	= 200	(ms)		
								
  Ra   = 0.1		
 	Rb   = 0.1		
 
 	celsius		(degC)
 } 
 
 
 ASSIGNED {
 	ik 		(mA/cm2)
	 ica (mA/cm2)
 	icaL (mA/cm2)
 	ek		(mV)
 	ninf
 	ntau 		(ms)	
  drive_channel1	(mM/ms)
  drive_channel2	(mM/ms)
 }
  
 
 STATE { 
 n 
 ca (mM)
	caL (mM)
}
 
 INITIAL { 
	ca=cainf
 caL=0
 cai=cainf
 rates(cai)
 	n = ninf
 }
 
 BREAKPOINT {
         SOLVE states METHOD cnexp
 	ik =  g *n* (v - ek)
 } 
 

DERIVATIVE states {  
 	
 	
 	
 	
 	
 	
 	

          rates(cai)    
         n' = (ninf-n)/ntau
}
PROCEDURE rates(cai(mM)) {  LOCAL a,b
							UNITSOFF
         a = Ra * (1e3*(cai  -cainf))^caix		
         b = Rb
         ntau = 1/(a+b)
        	ninf = a*ntau
					UNITSON
 }