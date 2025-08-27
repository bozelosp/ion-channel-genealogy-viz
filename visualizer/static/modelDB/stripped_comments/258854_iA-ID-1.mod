UNITS {

    (mA)    = (milliamp)
    (mV)    = (millivolt)
    (S)     = (siemens)
}





NEURON  {
    
        

    SUFFIX iA

    USEION k READ ek WRITE ik

        
        

    RANGE gkbar, gk, i

        

    GLOBAL minf, mtau, hinf

        

    THREADSAFE 

}





PARAMETER   {

    gkbar  = .0025      (S/cm2)   <0,1e9>
    ek     = -80        (mV)
    mtau    = 0.1       (ms)
}








ASSIGNED {
    
    v            (mV)
    celsius      (degC)

    gk           (S/cm2)
    ik           (mA/cm2)
    i            (mA/cm2)

    minf
    hinf
    htau         (ms)
}






STATE { 
    
    m h

}







? currents
BREAKPOINT  {

        

    SOLVE states METHOD cnexp

        

    gk = gkbar*m*h
    i  = gk*(v - ek) 
    ik = i

}





INITIAL     {

    rates(v)
    m = minf
    h = hinf

}








? states
DERIVATIVE states   {
    
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau

}





? rates
PROCEDURE rates( v(mV) ) {  

    LOCAL alpha, beta
    TABLE minf, hinf, htau FROM -100 TO 100 WITH 200   

        
        

        UNITSOFF

        
    
    alpha = 0.01 * vtrap( -(v+21.3) , 35.0 )
    beta  = 0.01 * vtrap(  (v+21.3) , 35.0 )
    minf  = alpha/(alpha+beta)

        
    
    alpha = -0.01 * vtrap( (v+58.0) , 8.2 )
    beta  = -0.01 * vtrap(-(v+58.0) , 8.2 )
    hinf  = alpha/(alpha+beta)

    htau = htauv(v)

        UNITSON

}



FUNCTION htauv( v (mV) ) {  

        UNITSOFF

        
    
    if ( v > - 20.0 ) {

            htauv = 5. + 2.6*( v+20. )/10.

    }else{

            htauv = 5.
    
    }

        UNITSON
}

FUNCTION vtrap(x,y) {  

        if (fabs(x/y) < 1e-6) {

                vtrap = y*(1.0 - x/y/2.0)

        }else{
        
                vtrap = x/(exp(x/y) - 1.0)
        }
}