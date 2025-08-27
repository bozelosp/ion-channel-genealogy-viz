UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
        (S) = (siemens)
        (mS) = (millisiemens)
}

NEURON {
    THREADSAFE
    

        SUFFIX KD2
        USEION k READ ek WRITE ik
        RANGE gmax, g
}

PARAMETER {

    gmax= 0.001 (S/cm2)
    np = 3		(1)

	kan=10.5	(mV)	
	Aan=3.0		(/ms)
	dan=-17		(mV)	
	kbn=30		(mV)	
	Abn=0.12	(/ms)
	dbn=-50		(mV)
	
	kal=35		(mV)	
	Aal=0.00025	(/ms) 
	dal=-85		(mV)	
	kbl=4.8		(mV)	
	Abl=0.0033	(/ms) 
	dbl=-51		(mV)	
}

ASSIGNED { 
    v (mV)
    ek (mV)
    
    ik (mA/cm2)
    an (/ms)
    bn (/ms)
    al (/ms)
    bl (/ms)
	g (S/cm2)
}

STATE {
	n
	l
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g  = gmax*n^np*l
    ik  = g*(v-ek)
}

INITIAL {
    settables(v)
    n = an/(an+bn)
    l = al/(al+bl)
}

DERIVATIVE states {  
    settables(v)      
    n' = ((an*(1-n)) - (bn*n))
    l' = ((al*(1-l)) - (bl*l))
}


PROCEDURE settables(v (mV)) {
    TABLE an, bn, al, bl DEPEND dan, dbn, dal, dbl, Aan, Abn, Aal, Abl
          FROM -100 TO 50 WITH 750

    if (v == dan) {
    	an = Aan
	} else {
		an = -Aan/kan*(dan-v) / (1 - exp((dan-v)/kan))	
	}
	bn = Abn*exp( (dbn-v)/kbn)	
	
	al = Aal*exp( (dal-v)/kal)	
	bl = Abl/(exp( (dbl-v)/ kbl) + 1)	
}