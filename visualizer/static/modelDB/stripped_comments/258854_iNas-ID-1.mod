UNITS {
    (mA) = (milliamp)
    (mV) = (millivolt)
	(S) = (siemens)
    R       = (k-mole) (joule/degC)     
    FARADAY = (faraday) (coulombs)      
}
 
? interface
NEURON {
    SUFFIX iNas
    USEION na READ ena WRITE ina
    RANGE gnabar, gna, Frec, i
    GLOBAL minf, hinf, sinf, mtau, htau, stau
	THREADSAFE 
}
 
PARAMETER {
    gnabar = .12 (S/cm2)	<0,1e9>
    Frec = 0.
}
 
STATE {
    m h s
}
 
ASSIGNED {
    v       (mV)
    celsius (degC)
    ena     (mV)

	gna     (S/cm2)
    ina     (mA/cm2)
    i       (mA/cm2)
    minf 
    hinf 
    sinf
	mtau    (ms) 
    htau    (ms) 
    stau    (ms)
}
 
? currents
BREAKPOINT {
    SOLVE states METHOD cnexp
    gna = gnabar*m*m*m*h*s
	i   = gna*(v - ena)
    ina = i
}
 
 
INITIAL {
	rates(v)
	m = minf
	h = hinf
    s = sinf
}

? states
DERIVATIVE states {  
    rates(v)
    m' = (minf-m)/mtau
    h' = (hinf-h)/htau
    s' = (sinf-s)/stau
}
 

? rates
PROCEDURE rates(v(mV)) {  
                  

    LOCAL  alpha, beta, Q
    TABLE minf, mtau, hinf, htau, sinf, stau DEPEND celsius FROM -100 TO 100 WITH 200

        UNITSOFF

    

            
    alpha = 0.182*vtrap( -(v+45), 4.5 )
    beta = 0.124*vtrap( (v+45), 4.5 )	
    mtau = 0.8/(alpha+beta)
    minf = alpha/(alpha+beta)

            
    alpha = 0.08*vtrap(-(v+40),3.0)
    beta  = 0.0005*vtrap( (v+10.0), 5. )
	htau = 1/(alpha+beta)
    hinf = 1.0 / ( 1 + exp((v+58)/5.) )

            
    Q = FARADAY / ( R * ( 273.16 + celsius ) )
    sinf = ( 1. + Frec*exp( (v+58.)/2. ) )/( 1. + exp( (v+58.)/2. ) )
    stau = 0.00333 (ms) * exp( 0.0024 (/mV) * (v+60.0) * Q )/( 1.0 + exp( 0.0012 (/mV) * (v+60.0) * Q )  )
}
 
FUNCTION vtrap(x,y) {  
    if (fabs(x/y) < 1e-6) {
            vtrap = y*(1 - x/y/2)
    }else{
            vtrap = x/(exp(x/y) - 1)
    }
}
 
        UNITSON