UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
	
}
 
NEURON {
 	SUFFIX Nam
	USEION na READ ena WRITE ina
	RANGE gnabar, gna, minf, hinf, ah, Bh
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	
	
	gnabar	= 0.035 (mho/cm2) 
	phi	= 5 < 0, 1e9 > 
	Vam	= -28 
	Kam	= 1
	Vbm	= -53 
	Kbm	= 18
	Vah	= -51 
	Kah	= 20
	Vbh	= -21 
	Kbh	= 1
	       
}
 
STATE {
        m h
}
 
ASSIGNED {
        ena (mV)
        v  (mV)
	ina (mA/cm2)
	celsius		(degC)
 	minf
	hinf
	ah
	Bh
        gna
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
        gna = gnabar*m^3*h
        ina = gna*(v - ena)
  
}
 
UNITSOFF
 
INITIAL {
	rates(v)
	m = minf
	h= hinf
}

DERIVATIVE states {  
        rates(v)      
       
	h'=phi*(ah*(1-h)-Bh*h)

}
 
PROCEDURE rates(v) {  
                      
        LOCAL  am, Bm
        
        
	am = (-0.1*(v-Vam)/Kam/(exp(-0.1*(v-Vam)/Kam)-1))
        Bm = 4*exp(-(v-Vbm)/Kbm)
        ah = 0.07*exp(-(v-Vah)/Kah)
        Bh=  1/(1+exp(-0.1*(v-Vbh)/Kbh))
        minf = am/(am+Bm)
 	hinf = ah/(ah+Bh)
	m=minf      
}
 
UNITSON