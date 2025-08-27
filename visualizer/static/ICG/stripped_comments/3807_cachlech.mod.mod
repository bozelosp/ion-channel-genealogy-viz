UNITS {  
        (mA) = (milliamp)  
        (mV) = (millivolt)  
}  
   
NEURON {  
        SUFFIX cach  
        USEION ca READ eca WRITE ica  
        RANGE gcabar 
        GLOBAL minf, mexp  
}  
   
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}  
   
PARAMETER {  
        v (mV)  
        celsius = 20 (degC)  
        dt (ms)  
        gcabar = 2e-06 (mho/cm2)  
        eca = 125 (mV)  
}  
   
STATE {  
        m   
}  
   
ASSIGNED {  
        ica (mA/cm2)  
        minf mexp  
}  
   
BREAKPOINT {  
        SOLVE states  
        ica = gcabar*m*(v - eca)  
}  
   
UNITSOFF  
   
INITIAL {  
     rates(v)  
     m = minf  
}  
PROCEDURE states() {  
        rates(v)      
        m = m + mexp*(minf-m)  
}  
   
PROCEDURE rates(v) {
         
         
         
     LOCAL  q10, tinc, alpha, beta, sum  
     TABLE minf, mexp DEPEND dt, celsius FROM -100 TO 100 WITH 200 

        q10 = 3^((celsius - 20)/10)  
        tinc = -dt * q10  
                
        alpha = 1.5 * vtrap(-(v-20),5)  
        beta =  1.5 * exp(-(v+25)/10)  
        sum = alpha + beta  
        minf = alpha/sum  
        mexp = 1 - exp(tinc*sum)  
}  
   
FUNCTION vtrap(x,y) {  
        if (fabs(x/y) < 1e-6) {  
                vtrap = y*(1 - x/y/2)  
        }else{  
                vtrap = x/(exp(x/y) - 1)  
        }  
}  
   
UNITSON  
  