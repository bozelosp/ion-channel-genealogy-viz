TITLE Cancer LP periodic PD input conductance

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
  SUFFIX synpd
  NONSPECIFIC_CURRENT i
  RANGE gbar, e, T, thresh, scale, phi, finf, f, g, i, vpre  : RANGE variables can vary over the length of a
                                                       : segment, and from segment to segment
  }

 
: All var names used must appear in some variable-declaration block,
: even if they're mentioned in the NEURON block already

PARAMETER  : this is a variable-declaration block, for params settable
           : from the interface
  {
  gbar = 0 (uS/mm2)
  e = -80 (mV)
  T = 1000 (ms)
  thresh = -58  (mV)
  scale = 15 (mV)
  tauf = 50 (ms)
  }

ASSIGNED  : also a var-declaration block, where one declares 1)
          : variables assigned to outside the mod file (not 
          : including parameters, I guess) , and 2) variables that
          : appear on the LHS of assignments within the mod file.
  {
  v (mV)
  phi
  finf
  vpre (mV)
  g (uS/mm2)
  i (mA/cm2)
  }

STATE  : also a var-decl block, where one declares state vars local to
       : this particular mechanism
  {
  f
  }

BREAKPOINT
  {
  SOLVE state_change METHOD cnexp  
  g = gbar*f
  i = (1e-4)*(g*(v-e))
  }

DERIVATIVE state_change
  {
  phi=2*PI*t/T
  vpre = (-62.6(mV)) + (11.039(mV))*cos(1*phi+(6.054)) + (5.199(mV))*cos(2*phi+(5.184)) + (1.437(mV))*cos(3*phi+(10.824)) + (0.323(mV))*cos(4*phi+(14.913)) + (0.256(mV))*cos(5*phi+(19.364)) + (0.325(mV))*cos(6*phi+(18.787)) + (0.110(mV))*cos(7*phi+(23.995))
    : These are the major Fourier components of a digitized _AB_ trace
  finf=naka_rushton((vpre-thresh)/scale)
  f' = (finf-f)/tauf
  }

INITIAL 
  {
  f = 0
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
 


  