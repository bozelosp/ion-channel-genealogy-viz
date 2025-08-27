NEURON {
        SUFFIX kdyn
        USEION k READ ko,ik WRITE ko 
        RANGE ko, ra, KAF, dep, peak
}

UNITS {
        (mM) = (milli/liter)
        (mA) = (milliamp)
        F    = (faraday) (coul)
}

PARAMETER {
        tck    = 1000   (ms)           
        koinf = 2 (mM)        
	kiinf = 140     (mM)	  
        dep   = 70e-3 (micron)     
	KAF   = 0.143 ()		  
peak = 3.03 ()
}

ASSIGNED {
        ik     (mA/cm2)
        ra
}

INITIAL {
	ko=koinf
}

STATE { ko (mM) 
}

BREAKPOINT { 
        SOLVE states METHOD derivimplicit
if (ko > peak) { ko = peak}
  if (ko < 2) { ko = 2}
}

DERIVATIVE states {      
 
        ko'= (1e4*(KAF*ik))/(F*dep)     
}