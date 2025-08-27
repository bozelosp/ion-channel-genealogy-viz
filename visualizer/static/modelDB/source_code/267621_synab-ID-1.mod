TITLE Cancer LP periodic AB input conductance

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
  PI = (pi) (1)
  }
 
NEURON  : this block manages the interface with NEURON
  {
  SUFFIX synab
  NONSPECIFIC_CURRENT i
  RANGE gbar, e, T, thresh, scale, phi, f, g, i,vpre
  }

PARAMETER
  {
  gbar = 0 (uS/mm2)
  e = -70 (mV)
  T = 1000 (ms)
  thresh = -58  (mV)
  scale = 15 (mV)
  }

ASSIGNED
  {
  v (mV)
  phi
  f
  vpre (mV)
  g (uS/mm2)
  i (mA/cm2)
  }

STATE
  {
  }

BREAKPOINT
  {
  phi=2*PI*(t/T)
  vpre = (-62.6(mV)) + (11.039(mV))*cos(1*phi+(6.054)) + (5.199(mV))*cos(2*phi+(5.184)) + (1.437(mV))*cos(3*phi+(10.824)) + (0.323(mV))*cos(4*phi+(14.913)) + (0.256(mV))*cos(5*phi+(19.364)) + (0.325(mV))*cos(6*phi+(18.787)) + (0.110(mV))*cos(7*phi+(23.995))
  f=naka_rushton((vpre-thresh)/scale)
  g = gbar*f
  i = (1e-4)*(g*(v-e))
  }

INITIAL 
  {
  }

FUNCTION naka_rushton(x)
  {
  LOCAL y
  if (x>0) {
    y=x^2
    naka_rushton=y/(1+y)
    }
  else {
    naka_rushton=0
    }
  }
 


  