UNITS
  {
  (pure) = (1)
  (mA) = (milliamp)
  (mV) = (millivolt)
  (mol) = (1)
  (nmol) = (nanomol)
  (molar) = (mol/liter)
  (mM) = (millimolar)
  (uM) = (micromolar)
  (uC) = (microC)
  (S) = (mho)
  (uS) = (microS)
  (pmol) = (picomole)
  (pC) = (picocoulomb)
  (J) = (joule)
  (um) = (micrometer)
  }

NEURON  
  {
  SUFFIX kca
  USEION k READ ek WRITE ik
  USEION ca READ cao, cai
  RANGE gbar, ik, g  
                     
  }

PARAMETER
  {
  gbar = 0      (uS/mm2)
  }

ASSIGNED
  {
  
  v (mV)
  cai (mM)
  cao (mM)
  celsius (degC)
  ek (mV)
  
  minf 
  taum (ms)
  hinf
  
  g (uS/mm2)
  ik (mA/cm2)
  }

STATE
  {
  m
  h
  }

BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  g = gbar*m*h
  ik = (1e-4)*(g*(v-ek))
  }

DERIVATIVE state_change
  {
  rates(v,cai)
  m' = (minf-m)/taum
  h' = (hinf-h)/11.85(ms)
  }

INITIAL
  {
  rates(v,cai)
  m=minf
  h=hinf
  }

PROCEDURE rates(v(mV),cai(mM)) 
  {  
  if (cai<1e-9)
    {
    cai=1e-9(mM)
    }
  minf=1/(1+(cai/1.43e-3(mM))^(-5))*1/(1+exp(-(v+5.5(mV))/8.0(mV)))
  taum=499.0(ms)+(5.0(ms)-499.0(ms))/
           (1+exp(-(v-(-51.9(mV)-(2.7(mV))*log(cai/(1e-3(mM)))))/10.0(mV)))
  hinf=1/(1+(cai/7.2e-3(mM))^1.25)
  }