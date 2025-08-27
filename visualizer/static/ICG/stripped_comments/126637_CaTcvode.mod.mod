UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX CaTcvode
        USEION ca READ eca WRITE ica
        RANGE  gcabar, ica, gca
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        
        gcabar = .0005 (mho/cm2)
        
   	
	 
	  mon = 1
	  hon = 1
	  alphaexp = 1
	  betaexp = 1
}
 
STATE {
        m h
}
 
ASSIGNED {
        eca (mV)
        ica (mA/cm2)
        gca minf tau q10 alpha beta sum hinf
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gca = gcabar * m*h
	ica = gca* (v-eca)
}
 
UNITSOFF
 
INITIAL {
	
	m = minf
	h = hinf
}

 
DERIVATIVE state {  
                      
        
        q10 = 3^((celsius - 37)/10)
        
                

        alpha = 2.6/(1+exp((v+21)/(-8)))
        beta =  0.18/(1+exp((v+40)/4))
        sum = alpha + beta
        minf = alpha/sum
	  tau= 1/(q10 * sum)
        m' = mon * (minf-m)/tau      

                

        alpha = 0.0025/(1+ (alphaexp * exp((v+40)/8)))
        beta = 0.19/(1+ (betaexp * exp((v+50)/(-10))))
        sum = alpha + beta
        hinf = alpha/sum
	  tau = 1/(q10 * sum)
	  h' = hon * (hinf-h)/tau
}

 
UNITSON