UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

 
? interface
NEURON {
        SUFFIX fuscurr
        USEION na READ ena WRITE ina
        USEION k READ ek WRITE ik
	USEION h READ eh WRITE ih VALENCE 1
        USEION nap READ enap WRITE inap VALENCE 1
        USEION kir READ ekir WRITE ikir VALENCE 1
	NONSPECIFIC_CURRENT il

	RANGE gna, gk, minf, hinf, ninf		
	GLOBAL mtau, htau, ntau, gnabar, gkbar	
		
	RANGE gh, kh_m_inf, kh_n_inf, aih	
	GLOBAL kh_m_tau, kh_n_tau, ghbar	

	RANGE gnap, gkir, minf_p, ninf_ir		
	GLOBAL mtau_p, ntau_ir, gnapbar, gkirbar	

	GLOBAL gl, el					
}

 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v (mV)
        celsius (degC)
        dt (ms)
	ek (mV)
        ena (mV)
        gnabar = 0.08 (mho/cm2)	<0,1e9>
        gkbar = 0.02 (mho/cm2)	<0,1e9>
        ghbar = 0.00054 (mho/cm2) <0,1e9>
	eh = -43 (mV)
	gl = 0.000150 (mho/cm2)	<0,1e9>
        el = -51.32 (mV) 
	mtau = 0.05 (ms) <0.01,100>
	htau = 0.5 (ms) <0.1,100>
	ntau = 0.5 (ms) <0.1,100>
	ekir (mV)
        enap (mV)
        gnapbar = 0.0001 (mho/cm2) <0,1e9>
        gkirbar = 0.0005 (mho/cm2)	<0,1e9>
	mtau_p = 5.3 (ms) <0,100>
	ntau_ir = 0.5 (ms) <0,100>
}


STATE {
        m h n khm khn mp nir
}
 
ASSIGNED {
    gna (mho/cm2)
    gk (mho/cm2)
    gh (mho/cm2)
    gnap (mho/cm2)
    gkir (mho/cm2)

    ina (mA/cm2)
    ik (mA/cm2)
    ih (mA/cm2)
    il (mA/cm2)
    inap (mA/cm2)
    ikir (mA/cm2)

    minf hinf
    ninf 
    kh_m_inf kh_n_inf
    kh_m_tau kh_n_tau	
    aih
    minf_p 
    ninf_ir 
}


LOCAL mexp, hexp, nexp, kh_m_exp, kh_n_exp, mexp_p, nexp_ir
 
? currents
BREAKPOINT {
    SOLVE states METHOD cnexp

    	gna = gnabar*m*m*h
	ina = gna*(v - ena)

   	gk = gkbar*n*n
	ik = gk*(v - ek)

	aih = 0.5*khm+0.5*khn
	gh = ghbar*aih
 	ih = gh*(v - eh)
	
	il = gl*(v - el)

    	gnap = gnapbar*mp*mp*mp
	inap = gnap*(v - enap)

   	gkir = gkirbar*nir
	ikir = gkir*(v - ekir)
}
? currents

UNITSOFF 
 

INITIAL {
	rates(v)
	m = minf
	h = hinf
	n = ninf
	khm = kh_m_inf
	khn = kh_n_inf
	mp = minf_p
	nir = ninf_ir
}


? states
DERIVATIVE states {  
    rates(v)
    m' = (minf - m) / mtau
    h' = (hinf - h) / htau
    n' = (ninf - n) / ntau
    khm' = (kh_m_inf - khm) / kh_m_tau
    khn' = (kh_n_inf - khn) / kh_n_tau
    mp' = (minf_p - mp) / mtau_p
    nir' = (ninf_ir - nir) / ntau_ir
}
 
LOCAL q10

? rates
PROCEDURE rates(v(mV)) {  
	                      
	LOCAL  alpha, beta, sum
	TABLE minf, mtau, hinf, htau, ninf, ntau, kh_m_inf, kh_n_inf, kh_m_tau, kh_n_tau, minf_p, ninf_ir DEPEND celsius FROM -200 TO 100 WITH 400

UNITSOFF
	q10 = 3^((celsius - 32)/10)

	
		minf = na_m(v)


	
		hinf = na_h(v)

	
        ninf = kd_m(v)

	
		kh_m_inf = kh_m(v) 
		kh_n_inf = kh_n(v)
		kh_m_tau = kh_mt(v)
		kh_n_tau = kh_nt(v) 

	
	minf_p = nap_m(v)

	
        ninf_ir = kird_m(v)
}

 
FUNCTION na_m(x) { 
	na_m = 1/(1+exp(-(x+38)/3.0))	
}

FUNCTION na_h(x) { 
	na_h = 1/(1+exp((x+43)/3.0))	
}

FUNCTION kd_m(x) { 
	kd_m = 1/(1+exp(-(x+40)/3))		
}

FUNCTION kh_m(x) { 
	kh_m = 1/(1+exp((x+87)/8.9))
}

FUNCTION kh_n(x) { 
	kh_n = 1/(1+exp((x+87)/8.9)) 
}

FUNCTION kh_mt(v) { 
	kh_mt = 100 + exp((v+183.6)/30.48)
}

FUNCTION kh_nt(v) { 
	kh_nt = 700 + exp((v+188.6)/11.2)/(1+exp((v+105)/5.5))
}

FUNCTION nap_m(x) { 
	nap_m = 1/(1+exp(-(x+55)/2.8))	
}

FUNCTION kird_m(x) { 
	kird_m = 1/(1+exp((x+85.48)/12)) 
}

UNITSON