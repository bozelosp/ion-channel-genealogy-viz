UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX CaP2cvode
        USEION ca READ eca WRITE ica
        RANGE  gcabar, ica, gca
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        gcabar = .0045 (mho/cm2)
        
	 
	 
	  mon = 1
}
 
STATE {
        m
}
 
ASSIGNED {
        eca (mV)
        ica (mA/cm2)
        gca minf tau q10 alpha beta sum
}
 
BREAKPOINT {
	  SOLVE state METHOD cnexp
        gca = gcabar * m
        ica = gca* (v-eca)
}
 
UNITSOFF
 
INITIAL {
	m = minf
}

DERIVATIVE state {   

                

        
        q10 = 3^((celsius - 37)/10)

        alpha = 8.5/(1+exp((v-8)/(-12.5)))
        beta =  35/(1+exp((v+74)/14.5))
        sum = alpha + beta
        minf = alpha/sum
	  tau = 1/(sum * q10)
	  m' = mon * (minf-m)/tau
}

 
UNITSON