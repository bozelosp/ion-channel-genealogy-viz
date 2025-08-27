UNITS {

    (mA)    = (milliamp)
    (mV)    = (millivolt)
    (S)     = (siemens)
    (molar) = (1/liter)
    (mM)    = (millimolar)
    R       = (k-mole) (joule/degC)     
    FARADAY = (faraday) (coulombs)      
}





NEURON  {
    
        

    SUFFIX iC

    USEION k READ ek WRITE ik
    USEION ca READ cai

        
        

    RANGE gkbar, gk, i

        

    GLOBAL minf, mtau

        

    THREADSAFE 

}





PARAMETER   {

    gkbar = .09075   (S/cm2)   <0,1e9>
    ek     = -80        (mV)
}








ASSIGNED {
    
    v            (mV)
    celsius      (degC)

    gk           (S/cm2)
    ik           (mA/cm2)
    i            (mA/cm2)

    cai          (mM)

    minf
    mtau         (ms)
}






STATE { 
    
    m 

}







? currents
BREAKPOINT  {

        

    SOLVE states METHOD cnexp

        

    gk = gkbar*m

    i  = gk*(v - ek)
    ik = i 

}





INITIAL     {

    rates(v,cai)
    m = minf

}








? states
DERIVATIVE states   {
    
    rates(v,cai)
    m' = (minf-m)/mtau

}





? rates
PROCEDURE rates( v(mV) , cai (mM) ) {  

    LOCAL alpha, beta, Q

        UNITSOFF

        
        

    Q = FARADAY / ( R * ( 273.16 + celsius ) )

    alpha = ( 0.48 (/ms) ) / ( 1 + ( 0.18 (mM) / cai )*exp( -1.68*v*Q ) )
    beta  = ( 0.28 (/ms) ) / ( 1 + cai/( 0.011 (mM) * exp( -2*v*Q ) ) )

        
    
    mtau = 1/(alpha+beta)
    minf = alpha/mtau

        UNITSON

}