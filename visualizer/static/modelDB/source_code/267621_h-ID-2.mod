TITLE Cancer LP Ih channel

COMMENT
ENDCOMMENT

UNITS
  {
  (mA) = (milliamp)
  (mV) = (millivolt)
  (molar) = (1/liter)
  (mM) = (millimolar)
  (S) = (mho)
  (uS) = (microS)
  }
 
NEURON
  {
  SUFFIX h
  NONSPECIFIC_CURRENT i
  RANGE e, gbar, g, i 
  }


PARAMETER  : this is a variable-declaration block, for params settable
           : from the interface
  {
  gbar = 0 (uS/mm2)
  e = -25 (mV)
  }

ASSIGNED 
  {
  v (mV)
  minf
  taum (ms)
  g (uS/mm2)
  i (mA/cm2)
  }

STATE
  {
  m
  }
 
BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  g = gbar*m
  i = (1e-4)*(g*(v-e))
  }

DERIVATIVE state_change
  {
  rates(v)  : Calculate minf, taum
  m' = (minf-m)/taum
  }

INITIAL 
  {
  rates(v)  : Calculate minf, taum
  m = minf
  }

PROCEDURE rates(v(mV)) 
  {  
  TABLE minf, taum FROM -150 TO 150 WITH 301
  minf=1/(1+exp( (v+84.3(mV))/6.4(mV)))
  taum=1/(1/(46.9(ms)*(1+exp(-(v- 29.7(mV))/19.4397(mV))))+
          1/(46.9(ms)*(1+exp( (v+206.2(mV))/19.4397(mV)))))
  }
