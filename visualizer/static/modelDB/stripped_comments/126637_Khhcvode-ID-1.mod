UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX Khhcvode
        USEION k WRITE ik
        RANGE   gk,  gkbar, ik
        
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius = 37 (degC)      
        gkbar = .036 (mho/cm2)
        ekcvode = -85(mV)
        non =1
}
 
STATE {
         n
}
 
ASSIGNED {
        ik (mA/cm2)
        gk ninf tau q10 alpha beta sum 
}
 
BREAKPOINT {
        SOLVE state METHOD cnexp 
        gk  = gkbar*n*n*n*n

        ik = gk*(v - ekcvode)      
}
 
UNITSOFF
 
INITIAL {
	
	n = ninf
}

DERIVATIVE state {  
                      
        

        q10 = 3^((celsius - 37)/10)
        
                
        alpha = .01*vtrapcvode(-(v+55),10) 
        beta = .125*exp(-(v+65)/80)
        sum = alpha + beta
        ninf = alpha/sum
        tau= 1/(q10 * sum)
        n' = non * (ninf-n)/tau      
}

FUNCTION vtrapcvode(x,y) {  
        if (fabs(x/y) < 1e-6) {
                vtrapcvode = y*(1 - x/y/2)
        }else{	
                vtrapcvode = x/(exp(x/y) - 1)
        }
}
 
UNITSON