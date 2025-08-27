INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS pGnoise
  POINTER vcell
  RANGE active
  RANGE g_e, g_i, E_e, E_i, g_e0, g_i0, g_e1, g_i1
  RANGE std_e, std_i, tau_e, tau_i, D_e, D_i
  RANGE del, dur, seed
  RANGE g_eout, g_iout
  USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
  (nA) = (nanoamp)
  (mV) = (millivolt)
  (umho) = (micromho)
}

PARAMETER {
  active = 0		
  
  dt (ms)

  del = 0	  (ms)
  dur = 5000 	  (ms)
  seed = 0
  
  E_e= 0 	  (mV)  
  E_i= -75 	  (mV)  

  g_e0= 0.0121    (umho)
  g_i0= 0.0573    (umho)

  std_e= 0.0030   (umho)
  std_i= 0.0066   (umho)

  tau_e= 2.728    (ms)  
  tau_i= 10.49    (ms)  
}

ASSIGNED {
  vcell	(mV)		
  iz    (nA)        	
  ez	(mV)
  g_e   (umho)          
  g_i   (umho)          
  g_e1  (umho)          
  g_i1  (umho)          
  D_e   (umho umho /ms) 
  D_i   (umho umho /ms) 
  exp_e
  exp_i
  amp_e (umho)
  amp_i (umho)

  g_eout (umho)         
  g_iout (umho)      	
}

INITIAL {
  set_seed(seed)
  
  exp_e = 0
  exp_i = 0
  amp_e = 0
  exp_i = 0

  g_e1 = 0
  g_i1 = 0
  
  if(tau_e != 0) {
    D_e = 2 * std_e * std_e / tau_e
    exp_e = exp(-dt/tau_e)
    amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )
  }
  
  if(tau_i != 0) {
    D_i = 2 * std_i * std_i / tau_i
    exp_i = exp(-dt/tau_i)
    amp_i = std_i * sqrt( (1-exp(-2*dt/tau_i)) )
  }
}

BREAKPOINT {
  if (active)  {
    if ((t >= del) && (t < dur)) {
      SOLVE oup
      
      if(tau_e==0) { 
        g_e1 = std_e * normrand(0,1) 
      }
      
      if(tau_i==0) { 
        g_i1 = std_i * normrand(0,1) 
      }
      
      g_e = g_e0 + g_e1
      if(g_e < 0) { g_e = 0 }
      
      g_i = g_i0 + g_i1
      if(g_i < 0) { g_i = 0 }
      
      iz = g_e * (vcell - E_e) + g_i * (vcell - E_i)

      g_eout = 50000.0 * g_e
      g_iout = 50000.0 * g_i
      
    } else { 
      iz = 0 
      g_eout = 0
      g_iout = 0
    }
  }
  else { 
    iz = 0 
    g_eout = 0
    g_iout = 0
  }
}

PROCEDURE oup() {       
  if(tau_e!=0) { 
    g_e1 =  exp_e * g_e1 + amp_e * normrand(0,1) 
  }
  
  if(tau_i!=0) { 
    g_i1 =  exp_i * g_i1 + amp_i * normrand(0,1) 
  }
}

PROCEDURE Update() {
  set_seed(seed)
  
  if(tau_e != 0) {
    D_e = 2 * std_e * std_e / tau_e
    exp_e = exp(-dt/tau_e)
    amp_e = std_e * sqrt( (1-exp(-2*dt/tau_e)) )
  }
  
  if(tau_i != 0) {
    D_i = 2 * std_i * std_i / tau_i
    exp_i = exp(-dt/tau_i)
    amp_i = std_i * sqrt( (1-exp(-2*dt/tau_i)) )
  }
}