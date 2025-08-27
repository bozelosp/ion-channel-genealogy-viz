INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
  POINT_PROCESS PRESYN
  RANGE spk_internal, spk
  GLOBAL thresh
}

PARAMETER {
  thresh = 0 			
}

ASSIGNED { 
  spk                           
  spk_internal                  
  v
}

INCLUDE "presyn.inc"

INITIAL {
  spk = 0
  spk_internal = 0
}

BREAKPOINT {
  SOLVE pp
}

PROCEDURE pp() {
  if (v > thresh) {    
    if (spk_internal == 0) {  
      newspike()                
      spk_internal = 1
      spk = 1
    }
  } else { 
    spk_internal = 0            
    spk = 0
  }
}