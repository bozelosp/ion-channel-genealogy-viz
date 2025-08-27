UNITS {

    (mA)    = (milliamp)
    (mV)    = (millivolt)
    (S)     = (siemens)
    R       = (k-mole) (joule/degC)     
    FARADAY = (faraday) (coulombs)      
}





NEURON  {
    
        

    SUFFIX iKDR

    USEION k READ ek WRITE ik

        
        

    RANGE gkbar, gk, i

        

    GLOBAL minf

        

    THREADSAFE 

}





PARAMETER   {

    gkbar  = .0014      (S/cm2)   <0,1e9>
    ek     = -80        (mV)
    mtau   = 3.50      (ms)
}








ASSIGNED {
    
    v            (mV)
    celsius      (degC)

    gk           (S/cm2)
    ik           (mA/cm2)
    i            (mA/cm2)

    minf
}






STATE { 
    
    m

}







? currents
BREAKPOINT  {

        

    SOLVE states METHOD cnexp

        

    gk = gkbar*m*m
    i  = gk*(v - ek) 
    ik = i

}





INITIAL     {

    rates(v)
    m = minf

}








? states
DERIVATIVE states   {
    
    rates(v)
    m' = (minf-m)/mtau

}





? rates
PROCEDURE rates(v(mV)) {  

    TABLE minf FROM -100 TO 100 WITH 200   

        
        

        UNITSOFF

        
    
    minf = 1.0/( 1.0 + exp( -(v+42.0)/3.0 ) )

        UNITSON

}