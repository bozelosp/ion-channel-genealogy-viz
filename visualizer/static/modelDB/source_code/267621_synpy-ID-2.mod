TITLE Cancer LP periodic PY input conductance

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
  SUFFIX synpy
  NONSPECIFIC_CURRENT i
  RANGE gbar, e, T, thresh, scale, phi, f, g, i, vpre  : RANGE variables can vary over the length of a
                                                 : segment, and from segment to segment
  }

 
: All var names used must appear in some variable-declaration block,
: even if they're mentioned in the NEURON block already

PARAMETER  : this is a variable-declaration block, for params settable
           : from the interface
  {
  gbar = 0 (uS/mm2)
  e = -70 (mV)
  T = 1000 (ms)
  thresh = -53  (mV)
  scale = 3 (mV)
  }

ASSIGNED  : also a var-declaration block, where one declares 1)
          : variables assigned to outside the mod file (not 
          : including parameters, I guess) , and 2) variables that
          : appear on the LHS of assignments within the mod file.
  {
  v (mV)
  phi
  f
  vpre (mV)
  g (uS/mm2)
  i (mA/cm2)
  }

STATE  : also a var-decl block, where one declares state vars local to
       : this particular mechanism
  {
  }

BREAKPOINT
  {
  phi=2*PI*t/T
  :vpre = (-57.6(mV)) + (11.425(mV))*cos(1*phi+(1.482)) + (4.407(mV))*cos(2*phi+(0.880)) + (1.581(mV))*cos(3*phi+(0.585)) + (1.177(mV))*cos(4*phi+(0.849)) + (0.786(mV))*cos(5*phi+(0.483)) + (0.598(mV))*cos(6*phi+(-6.109)) + (0.305(mV))*cos(7*phi+(-6.154))
    : these are the major Fourier components of a digitized PY trace
  vpre = (-57.6(mV)) + (11.425(mV))*cos(1*phi+(1.293)) + (4.407(mV))*cos(2*phi+(0.503)) + (1.581(mV))*cos(3*phi+(0.019)) + (1.177(mV))*cos(4*phi+(0.095)) + (0.786(mV))*cos(5*phi+(-0.460)) + (0.598(mV))*cos(6*phi+(-7.240)) + (0.305(mV))*cos(7*phi+(-7.474))
    : PY firing shifted later in the cycle
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
 


  