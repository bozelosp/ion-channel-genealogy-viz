NEURON {       
   SUFFIX ca       
   USEION ca READ ica, cai WRITE cai 
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
   cabath = 1.8   (mM)       
   diam = 1 (um)       
   ica            (mA/cm2)       
}       
       
STATE {       
   cai 
}       
       
       
INITIAL {       
   VERBATIM       
   cai = 0.0001;
   ENDVERBATIM       
}       
       
BREAKPOINT {       
   SOLVE state METHOD derivimplicit
}       
       
DERIVATIVE state {       
   cai' = -ica * 4/(diam*FARADAY) * (1e4)  
}