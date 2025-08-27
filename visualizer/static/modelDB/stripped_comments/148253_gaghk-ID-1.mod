NEURON {
	POINT_PROCESS gaghk
	USEION cl READ cli, clo WRITE icl VALENCE -1
	USEION hco3 READ hco3i, hco3o WRITE ihco3 VALENCE -1
	RANGE icl, ihco3, i

	RANGE GABAdur, GABA, taugaba, GABAINIT
	RANGE f1, f2
	RANGE kon, koff, koff2, k34, k43, alfa1, beta1, alfa2, beta2, alfa3, beta3 
	RANGE a1o, a1c, b1o, b1c, a2o, a2c, b2o, b2c, a3o, a3c, b3o, b3c
	RANGE d1, r1, d2, r2, d3, r3, d4, r4
	RANGE C1, C2, C3, C4, C5, C6, C7, C8, C9, C1O
	RANGE O1, O2, O3, D1, D2, D3, D4
	RANGE Prel, Pcl, Phco3, Rnumber
	RANGE ecl, ehco3, egaba
	RANGE gcl, ghco3, grel

	
}

UNITS {	
        (mA)    = (milliamp)
        (nA)    = (nanoamp)
	(mV)    = (millivolt)
	(uS)  = (micromho)
	(mM)    = (milli/liter)
	(uM)    = (micro/liter)
	F = (faraday) (coulombs)
	R = (k-mole)  (joule/degC)
}

PARAMETER {






	Prel  	= 0.18		      
	Rnumber	= 50		      

	Pcl    	= 8e-14 (cm3/s)	      

	


	celsius = 37    (degC)
	GABAdur	= 1	(ms)	      
	taugaba = 0.2   (ms)


	kon	= .003   (/uM /ms)	
	koff	= .150	 (/ms)		
	koff2   = .300   (/ms)		
	alfa1	= 1.111	 (/ms)		
	beta1	= .200	 (/ms)		
	alfa2	= .142   (/ms)		
	beta2	= 2.500  (/ms)		
	alfa3   = .150	 (/ms)		
	beta3   = .076   (/ms)		
	d1      = .013   (/ms)		
	r1	= .00013 (/ms)		
	d2	= .750   (/ms)		
	r2	= .015   (/ms)		
	d3	= .008   (/ms)		
	r3	= .00081 (/ms)		
	d4      = .00075 (/ms)		
	r4	= .00049 (/ms)		
	k34	= .710   (/ms)		
	k43	= .058   (/ms)		
	a1o	= .180	 (/ms)
	a2o	= .180	 (/ms)
	a1c	= 5.100  (/ms)
	a2c	= 5.100  (/ms) 
	b1o	= .070   (/ms)
	b2o	= .070   (/ms)
	b1c	= .630   (/ms)
	b2c	= .630   (/ms)
	a3o	= .090   (/ms)
	b3o	= .035   (/ms)
	a3c	= 5.100  (/ms)
	b3c	= .630   (/ms)
}


ASSIGNED {
	v		(mV)		

	cli	        (mM)
	clo		(mM)
	icl             (nA)		
	ecl		(mV)		

	hco3i		(mM)
	hco3o		(mM)
	ihco3           (nA)		
	ehco3		(mV)		

	egaba		(mV)		

	i               (nA)    	
	
	GABA		(mM)		
	GABAINIT	(mM) 		
	time0		(ms)
			
        f1		(/ms)    	
	f2		(/ms)   	
	
	
	Phco3  		(cm3/s)		
	gcl		(uS)		
	ghco3		(uS)		
	grel				

}

STATE {
        

	C1		
	C2		
	C3		
	C4		
	C5
	C6
	C7
	C8
	C9
	C1O
	O1		
	O2		
	O3		
	D1		
	D2		
	D3		
	D4		
}

INITIAL { 
	GABA		= 0
	GABAINIT	= 0
	C1		= 1
	C2		= 0
	C3		= 0
	C4		= 0
	C5		= 0
	C6		= 0
	C7		= 0
	C8		= 0
	C9		= 0
	C1O		= 0
	O1		= 0
	O2		= 0
	O3		= 0
	D1		= 0
	D2		= 0
	D3		= 0
	D4		= 0
	icl   = 0
	ihco3 = 0
	i     = 0

	Phco3 = Prel * Pcl	
}

BREAKPOINT {
	
	if (GABAINIT > 0) 
          
          { GABA = GABAINIT * exp(-(t-time0)/taugaba) }
	
	else {GABA = GABAINIT}
	
	SOLVE kstates METHOD sparse

	icl   = (1e+06)*(O1+O2+O3) * Pcl * Rnumber * ghk(v, cli, clo, -1)          
	ihco3 = (1e+06)*(O1+O2+O3) * Phco3 * Rnumber * ghk(v, hco3i, hco3o, -1) 
	i = icl + ihco3
	
	egaba = ghkvoltage(cli, clo, hco3i, hco3o)
	gcl = (1e+06)*(O1+O2+O3) * Pcl * Rnumber * conduct(v, cli, clo)		
	ghco3 = (1e+06)*(O1+O2+O3) * Phco3 * Rnumber *conduct(v, hco3i, hco3o)
	if (gcl>0) {grel = ghco3/gcl} else {grel = 0}
}

KINETIC kstates {
	f1    = 2 * kon * (1e3) * GABA
	f2    = kon * (1e3) * GABA

	~ C1 <-> C2	(f1,koff)
	~ C2 <-> C3	(f2,koff2)
	~ C2 <-> O1	(beta1,alfa1)
	~ C3 <-> O2	(beta2,alfa2)
	~ C2 <-> D1     (d1,r1)
	~ C3 <-> D2     (d2,r2)

	
	~ C3 <-> C4     (k34,k43)
	~ C4 <-> O3     (beta3,alfa3)
	~ C4 <-> D3     (d3,r3)
	~ D3 <-> D4     (d4,r4)

	~ O1 <-> C5     (a1c,a1o)
	~ O1 <-> C6     (b1c,b1o)
	~ O2 <-> C7     (a2c,a2o)
	~ O2 <-> C8     (b2c,b2o)
	~ O3 <-> C9     (a3c,a3o)
	~ O3 <-> C1O    (b3c,b3o)

	CONSERVE C1+C2+C3+C4+C5+C6+C7+C8+C9+C1O+O1+O2+O3+D1+D2+D3+D4 = 1
}


NET_RECEIVE(weight(mM), nspike) {

	
	if (flag == 0) {
		nspike = nspike + 1
		time0 = t
		GABAINIT  = weight
		
		net_send(GABAdur, nspike)
	}
	
	if (flag == nspike) {
		
		GABAINIT = 0
	}
}

FUNCTION ghk(v(mV), ci(mM), co(mM), z)  (millicoul/cm3) { 
	LOCAL e, w
        w = v * (.001) * z*F / (R*(celsius+273.15))
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else 
	
          { e = 1-w/2 }
        ghk = - (.001) * z* F * (co-ci*exp(w)) * e
}

FUNCTION ghkvoltage(c1i(mM), c1o(mM), c2i(mM), c2o(mM))  (mV) {

	ghkvoltage = - (1000)*(celsius + 273.15)*R/F*log((c1o + Prel*c2o)/(c1i + Prel*c2i))
}

FUNCTION conduct(v(mV), ci(mM), co(mM))  (millicoul/cm3/mV) { 
	LOCAL w
        w = v * (.001) *F / (R*(celsius+273.15))
	conduct = (0.001)*(.001)*F^2/(R*(celsius+273.15))*(ci-(co+ci)*exp(w)+(ci-co)*w*exp(w)+co*(exp(w)^2))/((1-exp(w))^2)
}