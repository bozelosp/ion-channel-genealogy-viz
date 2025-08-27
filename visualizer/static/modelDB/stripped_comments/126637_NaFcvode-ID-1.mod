UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX NaFcvode
	USEION na WRITE ina
        RANGE  gnabar, gna
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        
        gnabar	= 7.5 (mho/cm2)
        enacvode	= 45 (mV)
	mon = 1
	hon = 1
}
 
STATE {
        m h
}
 
ASSIGNED {
        ina (mA/cm2)
        gna minf hinf tau q10 alpha beta sum
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gna = gnabar *m*m* m*h 
	ina = gna* (v-enacvode)
}
 
UNITSOFF
 
INITIAL {
	
	m = minf
	h = hinf
}


DERIVATIVE state {  
                      
        
        q10 = 3^((celsius - 37)/10)
        
                

        alpha = 35/exp((v+5)/(-10))
        beta =  7/exp((v+65)/20)
        sum = alpha + beta
        minf = alpha/sum
	  tau= 1/(q10 * sum)
        m' = mon * (minf-m)/tau      

                

        alpha = 0.225/(1+exp((v+80)/10))
        beta = 7.5/exp((v-3)/(-18))
        sum = alpha + beta
        hinf = alpha/sum
    	  tau = 1/(q10 * sum)
        h' = hon * (hinf-h)/tau      
}

 
UNITSON