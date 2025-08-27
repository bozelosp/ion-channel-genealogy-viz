UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX Na
    USEION na READ nai,ena WRITE ina
    RANGE gna, m, h
    GLOBAL rest,rate_k
}

PARAMETER {
    v   (mV)
    dt  (ms)
    gna   = 1e-7    (mho/cm2)
    rest  = -60.0   (mV)
    ena             (mV)
    nai
    celsius	
}

STATE {
    m h  
}

ASSIGNED { 
    ina     (mA/cm2)
    alpham  (/ms)
    betam   (/ms)
    alphah  (/ms)
    betah   (/ms)
    rate_k
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    ina  = (gna)*m*m*h*(v-ena)
}

UNITSOFF

INITIAL {
    settables(v)
    m = alpham/(alpham+betam)
    h = alphah/(alphah+betah)
    rate_k = 1.78               
}

DERIVATIVE states {
    settables(v)      
    m' = alpham * (1-m) - betam * m
    h' = alphah * (1-h) - betah * h
}

PROCEDURE settables(v) {

    LOCAL vadj
    TABLE alpham, betam, alphah, betah DEPEND rest,celsius FROM -100 TO 100 WITH 400

    vadj = v - rest

    
    alpham = 0.32*vtrap((13.1-vadj),4)
    betam =  0.28*vtrap((vadj-40.1),5)  

    
    alphah = 0.128*exp((17-vadj)/18)
    betah = 4/(exp((40-vadj)/5)+1)

}

FUNCTION vtrap(x,y) {  
    if (fabs(x/y) < 1e-6) {
        vtrap = y*(1 - x/y/2)
    }else{
        vtrap = x/(exp(x/y) - 1)
    }
}
 
UNITSON