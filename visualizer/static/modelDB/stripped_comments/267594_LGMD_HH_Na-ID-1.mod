UNITS {
    (mV) = (millivolt)
    (mA) = (milliamp)
}

NEURON {
    THREADSAFE
    SUFFIX HH_Na
    USEION na READ ena WRITE ina
    RANGE gmax,g
}

PARAMETER {
    gmax= 0.045 (mho/cm2)

	kam=0.25
	Aam=3.0
	dam=-40
	kbm=0.05
	Abm=12
	dbm=-61

	kah=0.11
	Aah=2.25
	dah=-65
	kbh=0.1
	Abh=9.0
	dbh=-30
}

ASSIGNED { 
    v (mV)
    ena (mV)
    ina (mA/cm2)
    am (/ms)
    bm (/ms)
    ah (/ms)
    bh (/ms)
    g
}

STATE {
    m h
}

BREAKPOINT {
    SOLVE states METHOD cnexp
    g  = gmax*m^3*h
    ina  = g*(v-ena)
}

INITIAL {
    settables(v)
    m = am/(am+bm)
    h = ah/(ah+bh)
}

DERIVATIVE states {  
    settables(v)      
    m' = ((am*(1-m)) - (bm*m))
    h' = ((ah*(1-h)) - (bh*h))
}

UNITSOFF

PROCEDURE settables(v (mV)) {
    TABLE am, bm, ah, bh
          FROM -100 TO 100 WITH 2000

    if (v > (dbm-.1) && v < (dbm+0.1)) {
    	am = Aam
	} else {
		am = kam*-Aam*(dbm-v) / (1 - exp((dbm-v)*kam))
	}
	bm = Abm*exp( (dam-v)*kbm)
	
	ah = Aah*exp( (dah-v)*kah)
	bh = Abh/(exp( (dbh-v)* kbh) + 1)

}

UNITSON