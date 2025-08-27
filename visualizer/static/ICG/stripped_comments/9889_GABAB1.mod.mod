NEURON {  POINT_PROCESS GABAB1 }
















PARAMETER {
  Cmax    = 1     (mM)            
  Cdur    = 150   (ms)            
  Alpha   = 0.01  (/ms mM)        
  Beta    = 0.005 (/ms)           
  Erev    = -95   (mV)            
  Prethresh = 0                   
  Deadtime = 1    (ms)            
  gmax            (umho)          
}
INCLUDE "synq.inc"