UNITS {

    (mA)    = (milliamp)
    (mV)    = (millivolt)
    (S)     = (siemens)
    R       = (k-mole) (joule/degC)     
    FARADAY = (faraday) (coulombs)      
}





NEURON  {
    
        

    SUFFIX iM

    USEION k READ ek WRITE ik

        
        

    RANGE gkbar, gk, i, tau

        

    GLOBAL minf, mtau, Q, taua, taub

        

    THREADSAFE 

}





PARAMETER   {

    gkbar  = .06       (S/cm2)   <0,1e9>
    ek     = -80        (mV)
}








ASSIGNED {
    
    v            (mV)
    celsius      (degC)
    Q

    gk           (S/cm2)
    ik           (mA/cm2)
    i            (mA/cm2)

    minf
    mtau         (ms)
    taua         (ms)
    taub         (ms)
    tau          
}






STATE { 
    
    m 

}







? currents
BREAKPOINT  {

        

    SOLVE states METHOD cnexp

        

    gk = gkbar*m*Q*1e-4

    i  = gk*(v - ek) 
    ik = i

    tau = mtau

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

    LOCAL alpha, beta
    TABLE minf, mtau DEPEND celsius, Q FROM -100 TO 100 WITH 200   


        
        

        UNITSOFF

        

    Q = 2.3^((celsius - 23.)/10.)

    alpha = 1e-3 * vtrap(-(v+30.0), 9.0 )
    beta  = 1e-3 * vtrap( (v+30.0), 9.0 ) 

    mtau = 1.0 / ( Q * ( alpha + beta ) )
    minf = alpha / ( alpha + beta )
    
        UNITSON
}

FUNCTION vtrap(x,y) {  

        if (fabs(x/y) < 1e-6) {

                vtrap = y*(1 - x/y/2)

        }else{
        
                vtrap = x/(exp(x/y) - 1)
        }
}