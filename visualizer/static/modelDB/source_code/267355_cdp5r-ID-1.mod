COMMENT
  calcium accumulation into a volume of area*depth next to the
  membrane with a decay (time constant tau) to resting level
  given by the global calcium variable cai0_ca_ion
ENDCOMMENT

NEURON {
  SUFFIX cdp5r
USEION nca READ inca WRITE  ncai VALENCE 2
USEION lca READ ilca WRITE  lcai VALENCE 2
USEION tca READ itca WRITE tcai VALENCE 2
RANGE caiinf, caitot, ncai, lcai,tcai, ecatot, elca, enca, etca
GLOBAL catau
}

UNITS {
        (mV) = (millivolt)
  (molar) = (1/liter)
  (mM) = (milli/liter)
  (mA) = (milliamp)
  FARADAY = 96520 (coul)
  R = 8.3134  (joule/degC)
}

INDEPENDENT {t FROM 0 TO 100 WITH 100 (ms)}

PARAMETER {
        celsius = 6.3 (degC)
  depth = 200 (nm)  : assume volume = area*depth
  catau = 0.1 (ms)
  caiinf = 50.e-6 (mM)  : takes precedence over cai0_ca_ion
      : Do not forget to initialize in hoc if different
      : from this default.
  caotot = 2 (mM)
  ica (mA/cm2)
  inca (mA/cm2)
  ilca (mA/cm2)
  itca (mA/cm2)
  caitot= 50.e-6 (mM)
}

ASSIGNED {
  enca (mV)
  elca (mV)
  etca (mV)
  ecatot (mV)
}

STATE {
  ncai (mM)
  lcai (mM)
  tcai (mM)
}

INITIAL {
  VERBATIM
  ncai = _ion_ncai;
  lcai = _ion_lcai;
  tcai = _ion_tcai; 
  ENDVERBATIM
  ncai=caiinf/3
  lcai=caiinf/3
  tcai=caiinf/3
  caitot = caiinf  
  ecatot = ktf() * log(caotot/caiinf) 
  enca = ecatot
  elca = ecatot
  etca = ecatot
}


BREAKPOINT {
  SOLVE integrate METHOD derivimplicit
  caitot = ncai+lcai+tcai  
  ecatot = ktf() * log(caotot/caitot)  
  enca = ecatot
  elca = ecatot
  etca = ecatot
}

DERIVATIVE integrate {
ncai' = -(inca)/depth/FARADAY * (1e7) + (caiinf/3 - ncai)/catau
lcai' = -(ilca)/depth/FARADAY * (1e7) + (caiinf/3 - lcai)/catau
tcai' = -(itca)/depth/FARADAY * (1e7) + (caiinf/3 - tcai)/catau
}

FUNCTION ktf() (mV) {
  ktf = (1000)*R*(celsius +273.15)/(2*FARADAY)
} 