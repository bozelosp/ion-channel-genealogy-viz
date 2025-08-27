NEURON {
	SUFFIX bk
	USEION k READ ek WRITE ik
    RANGE minf, tm, ik
    RANGE gbar
    GLOBAL Vhalf, taumod
    GLOBAL vhm, vcm
    GLOBAL vhtm, atm, Ctm, tm0
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}

PARAMETER {

	ek		        (mV)

    gbar   = 0.05 	(mho/cm2)
    
    Vhalf   = -16   (mV)
    taumod  = 1

    vhm     = -16   (mV)
    vcm     = -8.5    (mV)

    vhtm    = -16 (mv)
    atm     = 10.0  (mV)
    Ctm     = 0.13  (mV)
    tm0     = 0.87 (ms)
}

STATE {
	m
}

ASSIGNED {
	v		        (mV)
	ik		(mA/cm2)
	minf
	tm		(ms)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
    
    
    ik = gbar * m * (v-ek)
}

DERIVATIVE states{
	rates(v)
	m' = (minf - m)/tm
}

INITIAL {
    setVhalf(Vhalf)
    setTauMod(taumod)
	rates(v)
	m = minf
}

PROCEDURE rates(v(mV)) {LOCAL q10
    q10 = 3^((celsius-34)/10)
    minf = 1/(1 + exp((v-vhm)/vcm))
    tm = (1/q10)*(tm0 + Ctm/(1 + exp((v-vhtm)/atm)))
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vhm = Vhalf
    vhtm = Vhalf
}

PROCEDURE setTauMod(taumod) {
    tm0 = 0.87 * taumod    
    Ctm = 0.13 * taumod
}