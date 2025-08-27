TITLE Cancer LP Ca current, drives KCa

COMMENT
ENDCOMMENT

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
  k = 1.381e-23 (J/K)
  q_e = 1.602e-19 (coulomb) 
  N_A = 6.022e23 (pure)
  z_ca = 2 (pure)
  }

NEURON  : this block manages the interface with NEURON
  {
  SUFFIX ca
  USEION ca READ cao, cai WRITE ica
  RANGE Pbar, h, frac, ica, ghknow
  }


PARAMETER  
  {
  Pbar = 0.0467 (nm/ms)
  }

ASSIGNED 
  {
  : things assigned to outside of this code
  v (mV)
  cao (mM)
  cai (mM)
  celsius (degC)
  ek (mV)
  : things used internally
  minf
  taum (ms)
  : things that get assigned to and then used in NEURON
  h
  frac
  ghknow (uC/cm2/nm)
  ica (mA/cm2)
  }

STATE
  {
  m
  }

BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  fh(cai)
  ghknow=ghk(v,cai,cao,z_ca)
  frac = m^3*h
  ica = Pbar*frac*ghknow
  }

DERIVATIVE state_change
  {
  rates(v)
  m' = (minf-m)/taum
  }

INITIAL
  {
  rates(v)
  m = minf
  }

PROCEDURE fh(cai(mM))
  {
  TABLE h FROM 0 TO 0.100 WITH 10001 
  h=1/(1+(cai/12.57e-3(mM)))
  }

PROCEDURE rates(v(mV)) 
  {  
  TABLE minf, taum FROM -150 TO 150 WITH 301
  minf=1/(1+exp(-(v+15.2(mV))/15.6(mV)))
  taum=1.8(ms)+(4.1(ms)-1.8(ms))/(1+exp(-(v+40.2(mV))/20.7(mV)))
  }

FUNCTION ghk(v(mV), ci(mM), co(mM), z) (uC/cm2/nm) 
  {
  LOCAL eta, eci, eco
  eta = (1e-3)*z*q_e*v/(k*(celsius+273.15))
  eco = co*efun(eta)
  eci = ci*efun(-eta)
  :high cao charge moves inward
  :negative potential charge moves inward
  ghk = (1e-7)*(N_A*q_e*z_ca*(eci - eco))
  }

FUNCTION efun(eta)
  {
  if (fabs(eta) < 1e-4)
    {
    efun = 1 - eta/2
    }
  else
    {
    efun = eta/(exp(eta) - 1)
    }
  }
