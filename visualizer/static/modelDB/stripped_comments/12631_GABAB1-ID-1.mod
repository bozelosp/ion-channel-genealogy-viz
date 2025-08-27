NEURON {  POINT_PROCESS GABAB1 }


















PARAMETER {
  Cdur    = 150   (ms)            
  Alpha   = 0.01  (/ms mM)        
  Beta    = 0.005 (/ms)           
  Erev    = -95   (mV)            
  Deadtime = 1    (ms)            
  GMAX     = 1    (umho)          
  DELAY = 0      (ms)            
}

INCLUDE "sns.inc"