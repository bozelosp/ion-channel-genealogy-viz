NEURON {
  POINT_PROCESS AMPA_S
  RANGE g, g_eff,g_specif,cellu_area
  RANGE Cdur, Alpha, Beta, Erev, Rinf, Rtau
  NONSPECIFIC_CURRENT i
}

UNITS {
  (nA)   = (nanoamp)
  (mV)   = (millivolt)
  (uS) 	= (microsiemens)
  (um)  = (microns)
}

PARAMETER {
  Cdur  = 1.0   (ms)  
  Alpha = 1.1   (/ms) 
  Beta  = 0.19  (/ms) 
  Erev  = 0     (mV)  
  g_specif = 1e-5	(uS/um2)		
  cellu_area= 1   (um2)
  g                 (uS)
}

ASSIGNED {
  v     (mV)   
  i     (nA)   
  g_eff (uS)   
  Rtau  (ms)   
  Rinf  
  synon 
}

STATE { Ron Roff }  




INITIAL {
  Rinf = Alpha / (Alpha + Beta)
  Rtau = 1 / (Alpha + Beta)
  synon = 0
  g=g_specif*cellu_area
}

BREAKPOINT {
  SOLVE release METHOD cnexp
  g_eff = g*(Ron + Roff)
  i = g_eff*(v - Erev)
}

DERIVATIVE release {
  Ron' = (synon*Rinf - Ron)/Rtau
  Roff' = -Beta*Roff
}

NET_RECEIVE(weight, on, r0, t0 (ms)) {
  
  
  if (flag == 0) {
    
    if (!on) {
      
      synon = synon + weight
      r0 = r0*exp(-Beta*(t - t0)) 
      
      
      
      
      
      Ron = Ron + r0
      Roff = Roff - r0
      t0 = t
      on = 1
      net_send(Cdur, 1)
    } else {
      
      net_move(t+Cdur)
    }
  }
  if (flag == 1) {
    
    
    synon = synon - weight
    
    
    r0 = weight*Rinf + (r0 - weight*Rinf)*exp(-(t - t0)/Rtau)
    
    
    
    Ron = Ron - r0
    Roff = Roff + r0
    t0 = t
    on = 0
  }
  
}