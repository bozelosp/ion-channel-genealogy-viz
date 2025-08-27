NEURON {POINT_PROCESS GABAB}
PARAMETER {
  Cdur    = 150   (ms)            
  Alpha   = 0.01  (/ms mM)        
  Beta    = 0.005 (/ms)           
  Erev    = -95   (mV)            
}
INCLUDE "netcon.inc"