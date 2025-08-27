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
    
        

    SUFFIX iAHPs

    USEION k READ ek WRITE ik
    USEION ca READ cai

        
        

    RANGE gkbar, gk, i

        

    GLOBAL minf, mtau, Q

        

    THREADSAFE 

}





PARAMETER   {

    gkbar   = .0005         (S/cm2)   <0,1e9>
    ek      = -80           (mV)
    caiNorm = 0.0001       (mM)
}








ASSIGNED {
    
    v            (mV)
    celsius      (degC)
    Q

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

        

    gk = gkbar*m*m*m
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
PROCEDURE rates( v(mV), cai(mM) ) {  

    LOCAL cac

        
        


        UNITSOFF

    cac = ( cai/( 0.00005 (mM) ) )^2
    Q   = 3^( (celsius-22)/10. )

        
    
    minf = cac/(1+cac)
    mtau = max( 0.5 , 0.003 (/ms) * ( 1+cac ) * Q )

        UNITSON

}





FUNCTION max( x, y ) {  

        UNITSOFF
    
    if ( x > y ) {

            max = x

    }else{

            max = y
    
    }

        UNITSON
}