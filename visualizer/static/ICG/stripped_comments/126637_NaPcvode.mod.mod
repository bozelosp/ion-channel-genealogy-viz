UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
 	SUFFIX NaPcvode
	USEION na READ ena WRITE ina
	RANGE gnabar, gna
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
  	v		(mV)
	celsius	= 37	(degC)
	
	gnabar	= 0.001 (mho/cm2)
	mon = 1
}
 
STATE {
        m 
}
 
ASSIGNED {
        ena (mV)
        ina (mA/cm2)
    	  minf
        gna
	  tau q10 alpha beta sum
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gna = gnabar*m*m*m
        ina = gna*(v - ena)
  
}
 
UNITSOFF
 
INITIAL {
	
	m = minf
}

 
DERIVATIVE state {  
                      
        
        q10 = 3^((celsius - 37)/10)
        
                

        alpha = 200/(1+exp((v-18)/(-16)))
        beta =  25/(1+exp((v+58)/8))
        sum = alpha + beta
        minf = alpha/sum
	  tau= 1/(q10 * sum)
        m' = mon * (minf-m)/tau      
}
 

 
UNITSON