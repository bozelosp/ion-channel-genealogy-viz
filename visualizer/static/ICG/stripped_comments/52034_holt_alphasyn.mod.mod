INDEPENDENT {t FROM 0 TO 10 WITH 1 (ms)}

NEURON {
  POINT_PROCESS holt_alphasyn
  RANGE i, stim, g, erev 
				
				
  NONSPECIFIC_CURRENT i
}

PARAMETER {
  stim = 0 (uS)			
				
  erev = 0 (mV)			
  v (mV)
  dt (ms)
}

ASSIGNED {
  tau (ms)			
  i (nA)			
  state_t1 (uS)			
  state_t2 (uS)			
				
  new_state
  scale_factor			
  A_coeff			
  B_coeff			
  g (uS)
}

INITIAL {
  IF (tau == 0.0) {
    tau = 1			
  }
  reinit()			
}





PROCEDURE set_tau(new_tau) {
  tau = new_tau			
  reinit()			
}

FUNCTION get_tau() {		
  get_tau = tau
}
































































PROCEDURE reinit() {
  state_t1 = 0			
  state_t2 = 0			
  A_coeff = 2 * exp(-dt/tau)	
  B_coeff = exp(-2 * dt/tau)
  scale_factor = (dt/tau)*exp(dt/tau+1) 
				
}




BREAKPOINT {

SOLVE dum
  g = new_state
  i = g*(v-erev)
}

PROCEDURE dum() {
    state_t2 = state_t2 - stim*scale_factor
    stim = 0
    new_state = A_coeff * state_t1 - B_coeff * state_t2
    state_t2 = state_t1		
    state_t1 = new_state

VERBATIM
return 0;
ENDVERBATIM
  }