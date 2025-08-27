INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS mGLUR
  RANGE G, C, lastrelease
  RANGE Cmax, Cdur, Deadtime, K1, K2, initmGluR
  RANGE rip3 
  RANGE ip3 
  RANGE scalef 
  RANGE degGluRate
  
}

UNITS {
  (mM) = (milli/liter)
  (uM)=  (micro/liter)
  (mA)    = (milliamp) 
}

PARAMETER {	
  initmGluR=0.3e-3 (mM)
  Cmax	= 1	(mM)		
  Deadtime = 1	(ms)		
  K1	= 0.28	(/ms uM)	
  K2	= 0.016 (/ms)		
  K_PLC = 5 (uM)		
  K_PIP2 = 160 (uM)		
  K_G=25 (uM)
  
  kfplc = 0.83(/ms)
  kbplc = 0.68 (/ms) 
  Vmax1 = 0.58 (/ms)
  
  D5f = 15 (/ms)
  D5b = 7.2 (/ms)
  D6f = 1.8(/ms)
  Vmax2 = 1.8 (/ms)
  Km2 = 0.6 (uM)
  
  G2f = 100 (/ms)
  G2b = 100 (/ms)
  
  D7f =9  (/ms)
  G9f = 0.75(/ms)  
  Cdur=2 (ms)			
  scalef = 1e10
  degGluRate = 1.0 (/ms) 
}

ASSIGNED {
  C		(mM)		
  lastrelease	(ms)		
  rip3
}

STATE {
  aG				
  aPLC_aG
  aPLC_PIP2
  Glu 
  degGlu 
  Glu_mGluR
  GG_mGluR
  ip3
  
  mGluR
  PLC
  PIP2
  G
  ip3i
}

INITIAL {
  Glu = 0
  degGlu = 0
  Glu_mGluR = 0
  GG_mGluR = 0
  aPLC_aG=0 
  aPLC_PIP2=0
  aG =0
  ip3= 0
  G=K_G-(aG+GG_mGluR+aPLC_aG+aPLC_PIP2)
  PLC=K_PLC-(aPLC_aG+aPLC_PIP2)
  PIP2=K_PIP2-(aPLC_PIP2+ip3)
  mGluR=initmGluR
  lastrelease = -1e8
  rip3 = 0
  
}

BREAKPOINT {
  
  SOLVE bindkin METHOD sparse
}

KINETIC bindkin { LOCAL a,b

  

  ~ Glu+mGluR <-> Glu_mGluR (K1, K2)
  ~ Glu <-> degGlu (degGluRate,0)

  ~ Glu_mGluR + G <-> GG_mGluR (D5f,D5b)
  ~ GG_mGluR <-> aG+mGluR (D6f,0)
  ~ aG <-> G (D7f,0)
  ~ aG+PLC <-> aPLC_aG (G2f, G2b)
  ~ aPLC_aG+PIP2 <-> aPLC_PIP2 (kfplc,kbplc)
  ~ aPLC_PIP2 <-> ip3 (Vmax1,0)
  rip3 = f_flux * scalef 
  
  
  
  
  
  
  
}

PROCEDURE evaluateC() {
  LOCAL q
  q = ((t - lastrelease) - Cdur)		
  if (q >= 0 && q <= Deadtime && C == Cmax) {	
    C = 0.
  }
}

NET_RECEIVE (weight)  { 
  Glu = Glu + weight 









}