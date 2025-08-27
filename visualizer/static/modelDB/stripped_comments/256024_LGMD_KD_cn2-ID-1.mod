UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(S) = (siemens)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

NEURON {
    THREADSAFE
    

        SUFFIX KD_cn2
        USEION k READ ek WRITE ik
        USEION cn READ cni VALENCE 1
        RANGE gmax, g, taun, taul
        GLOBAL tnmax, tlmax
}

PARAMETER {
	gmax=0.01 (mho/cm2)
	vhalfn=-52	(mV)
	vn2=-65		(mV)
	zn=6.0		(mV)
	tnmax=20	(ms)
	tnmin=1.5	(ms)
	tns=-7.5	(mV)
	np=2
	
	vhalfl=-67	(mV)
	zl=-3.6		(mV)
	tlmax=2000	(ms)
	tlmin=50	(ms)
	vl2=-75		(mV)
	tls=14		(mV)
	
	cnvm=25		(mV)
	lcp=1		(1)
	kD=3e-4		(mM)
	taucn=20	(ms)
}

STATE {
	n
	l
	ov	(mV)
	ovs	(mV)
}

ASSIGNED {
    v (mV)
    ek (mV)
    cni (mM)

	ik (mA/cm2)
	ninf
	linf
	vs	(mV)
	taul (ms)
	taun (ms)
	g (S/cm2)
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	l = 1/(1+exp((vhalfl-ov+ovs)/zl))
	g = gmax*n^np*l
	ik = g*(v-ek)
}

INITIAL {
	lf(v)
	rates(v)
	n = ninf
	l = 1/(1+exp((vhalfl-v)/zl)*exp(vs/zl))
	ov = v
	ovs = vs
}

FUNCTION alpn(v(mV)) {
  alpn = exp((vhalfn-v)/zn)
}

FUNCTION betn(v(mV)) {
  betn = exp((vn2-v)/tns) 
}

DERIVATIVE states {
	lci(cni)
	rates(v)
	n' = (ninf - n)/taun
	
	
	
	ov' = (v-ov)/taul
	ovs' = (vs-ovs)/taucn
}

PROCEDURE lf(v(mV)) {
	LOCAL dvdt, dsdt
	lci(cni)
	
	linf = 1/(1+exp((vhalfl-ov+ovs)/zl))
}

PROCEDURE rates(v (mV)) { 
	LOCAL a
	TABLE ninf, taun, taul DEPEND vhalfn, tlmax, tnmax, tnmin
          FROM -100 TO 50 WITH 600
    
	a = alpn(v)
	ninf = 1/(1 + a)
	taun = 4*(tnmax-tnmin)/(1+betn(v))*ninf+tnmin
	taul = 2*tlmax/( exp((v-vl2)/tls) + exp((vl2-v)/tls) ) + tlmin
}


PROCEDURE lci(cni (mM)) { 
	TABLE vs DEPEND lcp, kD, cnvm
          FROM 0 TO 0.01 WITH 500
    
	vs = cnvm-cnvm/(1+(cni/kD)^lcp)
    
}