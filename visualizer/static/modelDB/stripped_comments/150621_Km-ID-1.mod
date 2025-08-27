UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
 	SUFFIX Km
	USEION k WRITE ik
	RANGE gkmbar, gkm, ninf, an, Bn
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	
	ek	= -90	(mV)
	gkmbar	= 0.006 (mho/cm2) 
	phi	= 5  < 0, 1e9 > 
      Van	= -27 
	Kan	= 1
	Vbn	= -37 
	Kbn	= 80
}
 
STATE {
        n
}
 
ASSIGNED {
        v  (mV)
	ik (mA/cm2)
	celsius		(degC)
 	ninf
	an
	Bn
        gkm
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gkm = gkmbar*n^4
        ik = gkm*(v - ek)
  
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	n= ninf
}

DERIVATIVE states {  
        rates(v)      
       
	n'=phi*(an*(1-n)-Bn*n)

}
 
PROCEDURE rates(v) {  
                      
        
        an = (-0.01*(v-Van)/Kan/(exp(-0.1*(v-Van)/Kan)-1))
        Bn = 0.125*exp(-(v-Vbn)/Kbn)
        ninf = an/(an+Bn)
	      
}
 
UNITSON