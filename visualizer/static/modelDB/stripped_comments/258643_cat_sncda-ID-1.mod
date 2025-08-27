INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
	SUFFIX cat
	USEION ca READ eca,cai,cao WRITE ica
	RANGE m, h, hs, gcat, gbar, ica
    EXTERNAL apc_metap, fpc_metap
	GLOBAL vhm, vhh, vhhs, km, kh, khs
    GLOBAL tm0, th0, ths0, vhtm, vhth, atm, ath, Ctm, Cth
	RANGE minf, hinf, hsinf, mtau, htau, hstau
	GLOBAL q10, tadj
    GLOBAL Vhalf, taumod 
    
    
    
}

PARAMETER {
	gbar = 1   	(cm/s) 
	Vhalf = -40	(mV)		
    taumod = 1
    cai         (mM)
    cao         (mM)
    
    
    
    
    
    
    
    

    vhm  = -54.5 (mV)    
    km   = 5   (mV)    

	vhh  = -64.5	(mV)	
    kh   = -1.6 (mV)        
    
    vhhs = -64.5	(mV)	
    khs  = -1.6    (mV)        

    tm0  = 3.2
    th0  = 76
    ths0 = 600

    vhtm  = -40 (mV)
    vhth  = -46 (mV)

    atm   = 4.6
    ath   =8.85

    Ctm   = 19
    Cth   = 43

	temp = 33	(degC)		
	q10  = 3.0			

	v 		(mV)
	dt		(ms)
	celsius		(degC)

    tadj = 1

}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
    F = 9.6485e4    (coul)
    R = 8.3145      (joule/degC)
	(S) = (siemens)
	(mM) = (milli/liter)
	(pS) = (picosiemens)
	(um) = (micron)
} 

ASSIGNED {
	ica 		(mA/cm2)
	gca		(cm/s) 
	eca		(mV)	
	minf 		hinf        hsinf
	mtau (ms)	htau (ms)   hstau (ms)
    T   (degC)
    E   (volts)
    z
}
 

STATE { m h hs }

INITIAL { 
    
    setVhalf(Vhalf)
    setTauMod(taumod)
	m = minf_cat(v)
	h = hinf_cat(v)
    hs = hsinf_cat(v)
    tadj = tadj_ca_t()
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gca = gbar*m*m*m*((0.6*h)+(0.4*hs))
	    ica = gca * ghk(v,cai,cao) 
} 

FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (degC)) (mV) {
        KTF = ((25.26 (mV) /293.15 (degC) )*(celsius + 273.15 (degC) ))
}

FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

DERIVATIVE states {   
        m' = -(m-minf_cat(v))/(tadj*taum_cat(v))
        h' = -(h-hinf_cat(v))/(tadj*tauh_cat(v))
        hs' = -(hs-hsinf_cat(v))/(tadj*tauhs_cat(v))
}

FUNCTION boltz(x,y,z) {
		boltz = 1/(1+exp(-(x-y)/z))
}

FUNCTION minf_cat(v (mV)) (1) {
        minf_cat = boltz(v,vhm,km) 
}

FUNCTION taum_cat(v (mV)) (1/ms) {
        taum_cat = tm0 + Ctm/(1+exp((v-vhtm)/atm))
}

FUNCTION hinf_cat(v (mV)) (1) {
        hinf_cat = (boltz(v,vhh,kh))
}

FUNCTION tauh_cat(v (mV)) (1/ms) {
        tauh_cat = th0 + Cth/(1+exp((v-vhth)/ath))
}

FUNCTION hsinf_cat(v (mV)) (1) {
        hsinf_cat = (boltz(v,vhhs,khs))
}

FUNCTION tauhs_cat(v (mV)) (1/ms) {
        tauhs_cat = ths0
}

FUNCTION tadj_ca_t() {
        tadj_ca_t =  1/(q10^((celsius - temp)/10))
}

PROCEDURE setVhalf(Vhalf(mV)) {
    vhm = Vhalf
    vhh = Vhalf - 10
    vhhs = Vhalf - 10
    vhtm = Vhalf + 14.5
    vhth = Vhalf + 8.5
}

PROCEDURE setTauMod(taumod) {
    tm0 = 3.2 * taumod
    th0 = 76 * taumod
    ths0 = 600 * taumod
    Ctm = 19 * taumod
    Cth = 43 * taumod
}