NEURON {   
   SUFFIX na   
   USEION na READ ina, nai, nao WRITE nai, nao   
}   
   
INDEPENDENT {t FROM 0 TO 1 WITH 10 (ms)}   
   
UNITS {   
   (molar) = (1/liter)   
   (mV) = (millivolt)   
   (um) = (micron)   
   (mM) = (millimolar)   
   (mA) = (milliamp)   
   FARADAY = 96520 (coul)   
   R = 8.3134     (joule/degC)   
}   
   
PARAMETER {   
   celsius=20     (degC)   
   nabath = 115   (mM)   
   diam = 1 (um)   
   ina            (mA/cm2)   
}   
   
STATE {   
   nai START 10 (mM)   
   nao START 115  (mM)   
}   
   
   
INITIAL {   
   VERBATIM   
   nai = _ion_nai;   
   nao = _ion_nao;   
   ENDVERBATIM   
}   
   
BREAKPOINT {   
   SOLVE state METHOD derivimplicit
}   
   
DERIVATIVE state {   
   nai' = -ina * 4/(diam*FARADAY) * (1e4)
   nao'= 0
}   
      
   
