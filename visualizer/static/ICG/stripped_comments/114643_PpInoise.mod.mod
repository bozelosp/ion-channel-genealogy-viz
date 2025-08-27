INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS pInoise
  RANGE active
  RANGE I_0, std_I, tau_I, D_I
  RANGE del, dur, seed
  USEION z READ ez WRITE iz VALENCE 1
}

UNITS {
  (nA) = (nanoamp)
}

PARAMETER {
  active = 0		
  
  dt              (ms)

  del = 0	  (ms)
  dur = 5000 	  (ms)
  seed = 0
  
  I_0    = 0.1    (nA)  
  std_I  = 0.1    (nA)  
  tau_I  = 2.0    (ms)  
}

ASSIGNED {
  iz    (nA)        	
  ez	(mV)
  iz1   (nA)            
  D_I   (nA nA /ms)     
  exp_I
  amp_I (nA)
}

INITIAL {
  set_seed(seed)
  
  iz1 = 0
  
  if(tau_I != 0) {
    D_I = 2 * std_I * std_I / tau_I
    exp_I = exp(-dt/tau_I)
    amp_I = std_I * sqrt( (1-exp(-2*dt/tau_I)) )
  }
}

BREAKPOINT {
  if (active)  {
    if ((t >= del) && (t < dur)) {
      SOLVE oup
      iz = - I_0 - iz1
    } else { iz = 0 }
  }
  else { iz = 0 }   
}

PROCEDURE oup() {       
  if (tau_I!=0) { iz1 = exp_I * iz1 + amp_I * normrand(0,1) }
  if (tau_I==0) { iz1 = std_I * normrand(0,1) }
}

PROCEDURE Update() {
  set_seed(seed)
  
  if(tau_I != 0) {
    D_I = 2 * std_I * std_I / tau_I
    exp_I = exp(-dt/tau_I)
    amp_I = std_I * sqrt( (1-exp(-2*dt/tau_I)) )
  }
}