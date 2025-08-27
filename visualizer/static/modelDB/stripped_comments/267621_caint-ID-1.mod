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

NEURON
  {
  SUFFIX caint
  USEION ca READ ica WRITE cai
  GLOBAL Pbar1, vol1
  RANGE Pbar
  }


PARAMETER 
  {
  Pbar = 0.0467 (nm/ms)
  Pbar1 = 1.1675 (um3/s)
  vol1 = 6.49 (um3)
  }

ASSIGNED 
  {
  ica (mA/cm2)
  }

STATE 
  {
  cai (mM)
  }

BREAKPOINT
  {
  SOLVE state_change METHOD cnexp
  }

DERIVATIVE state_change
  {
  cai' = (0.020e-3(mM)-cai)/70.4(ms) +
         (1e4)*((-1/(z_ca*q_e*N_A*vol1))*Pbar1/Pbar*ica)
  }

INITIAL
  {
  cai = 0.57e-3(mM)  
  }