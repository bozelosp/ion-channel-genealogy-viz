UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX CaEcvode
        USEION ca READ cai, cao WRITE ica
        RANGE  gcabar, ica, gca
} 
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)
        mon = 1
	  hon = 1
        gcabar = .0005 (mho/cm2)
        ecacvode = 135 (mV)
   	  cai	= 0.40e-4 (mM)		
	  cao	= 2.4	(mM)

}
 
STATE {
        m h
}
 
ASSIGNED {
        ica (mA/cm2)
        gca minf tau q10 alpha beta sum hinf
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp
        gca = gcabar * m*h
	ica = gca* (v-ecacvode)
}
 
UNITSOFF
 
INITIAL {
	
	m = minf
	h = hinf
}


DERIVATIVE state {  
                      
        
        q10 = 3^((celsius - 37)/10)
        
                

        alpha = 2.6/(1+exp((v+7)/(-8)))
        beta =  0.18/(1+exp((v+26)/4))
        sum = alpha + beta
        minf = alpha/sum
	  tau= 4/(q10 *sum)
        m' = mon * (minf-m)/tau      

                

        alpha = 0.0025/(1+exp((v+32)/8))
        beta = 0.19/(1+exp((v+42)/(-10)))
        sum = alpha + beta
        hinf = alpha/sum
	  tau = 10/(q10 * sum)
	  h' = hon * (hinf-h)/tau
}

 
UNITSON