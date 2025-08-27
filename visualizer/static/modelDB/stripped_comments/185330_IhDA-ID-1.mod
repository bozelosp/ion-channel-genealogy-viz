UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	
}

PARAMETER {
	v 		(mV)
      ehd=-30  		(mV)        
	ghdbar=0.001	(mho/cm2) 
	vh=-95 (mV) 
	tc 
        vhalft=-112   	(mV)
        a0t=0.0016      	(/ms)
        zetat=2.2    	(1)
        gmt=.9   	(1)
	k=8
}

NEURON {
	SUFFIX hd
	NONSPECIFIC_CURRENT i
        RANGE ghdbar
        GLOBAL linf,taul, vh, tc
}

STATE {
        l
}

ASSIGNED {
	  i (mA/cm2)
        linf      
        taul
        ghd
}

INITIAL {
	rate(v)
	l=linf
}


BREAKPOINT {
	SOLVE states METHOD cnexp
	ghd = ghdbar*l
	i = ghd*(v-ehd)

}

DERIVATIVE states {     
        rate(v)
        l' =  (linf - l)/taul
}

FUNCTION alpt(v(mV)) {
  alpt = exp(0.0378*zetat*(v-vhalft)) 
}

FUNCTION bett(v(mV)) {
  bett = exp(0.0378*zetat*gmt*(v-vhalft)) 
}


PROCEDURE rate(v) { 
	LOCAL a
        linf = 1/(1 + exp((v-vh)/k))

        a = alpt(v)
        taul = bett(v)/(a0t*(1+a))
}