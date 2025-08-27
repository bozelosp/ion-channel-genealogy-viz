UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {
	v 		(mV)
	celsius 	(degC)
	gbar	= .003 	(mho/cm2)
    vhalfn	= -3.8   	(mV)
    a0n	= 0.02  (/ms)
    zetan	= -3    (1)
    gmn	= 0.5  	(1)
	nmax	= 2  	(1)
	q10	= 1	(1)
	FRT = 39 (coulombs/joule) 
}

NEURON {
	SUFFIX Kdr
	USEION k WRITE ik
    RANGE  gkdr,gbar,ik
	GLOBAL ninf,taun
}

STATE { n }

ASSIGNED {
	ik 	(mA/cm2)
    ninf	
    gkdr
    taun
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	gkdr 	= gbar * n
	ik 	= gkdr * ( v + 90.0 )
}

INITIAL {
	rates(v)
	n=ninf
}

FUNCTION alpn(v(mV)) {
  alpn = exp( 1.e-3 * zetan * ( v - vhalfn ) * FRT ) 
}

FUNCTION betn(v(mV)) {
  betn = exp( 1.e-3 * zetan * gmn * ( v - vhalfn ) * FRT ) 
}

DERIVATIVE states {     
        rates(v)
        n' = ( ninf - n ) / taun
}

PROCEDURE rates(v (mV)) { 
        LOCAL a,qt
        qt	= q10 ^ ( ( celsius - 24 ) / 10 )
        a 	= alpn(v)
        ninf 	= 1 / ( 1 + a )
        taun 	= betn(v) / ( qt * a0n * ( 1 + a ) )
	if (taun<nmax) { taun=nmax }
}