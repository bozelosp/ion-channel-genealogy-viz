NEURON{
	SUFFIX SK2
	USEION ca READ cai
	USEION k READ ek WRITE ik 
	USEION lca READ lcai VALENCE 0
	RANGE gkbar, gk, ik, acai
	GLOBAL diff, Q10, fac, invc1,invc2,invc3,invo1,invo2,diro1,diro2,dirc2,dirc3,dirc4
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
	gkbar = 0.038 (mho/cm2)
	Q10 = 3 (1)
	diff = 3 (1) 
	fac = 1 
	ocontrt = 0

	invc1 = 0.08  ( /ms)
	invc2 = 0.08  ( /ms)
	invc3 = 0.2 ( /ms)

	invo1 = 1.0     ( /ms)
	invo2 = 0.1     ( /ms)
	diro1 = 0.16 ( /ms)
	diro2 = 1.2    ( /ms)


	dirc2 = 200 ( /ms-mM )
	dirc3 = 160 ( /ms-mM )
	dirc4 = 80  ( /ms-mM )

}

ASSIGNED{ 
	v	(mV) 
	ek	(mV) 
	cai (mM)
	celsius		(degC)
	gk	(mho/cm2) 
	ik	(mA/cm2) 
	invc1_t  ( /ms)
	invc2_t  ( /ms)
	invc3_t  ( /ms)
	invo1_t  ( /ms)
	invo2_t  ( /ms)
	diro1_t  ( /ms)
	diro2_t  ( /ms)
	dirc2_t  ( /ms-mM)
	dirc3_t  ( /ms-mM)
	dirc4_t  ( /ms-mM)
	tcorr	 (1)

	dirc2_t_ca  ( /ms)
	dirc3_t_ca  ( /ms)
	dirc4_t_ca  ( /ms)
	
	lcai		(mM)
	acai     (mM)
	
	
	
} 

STATE {
	c1
	c2
	c3
	c4
	o1
	o2
}

BREAKPOINT{ 
	SOLVE kin METHOD sparse 
	gk = gkbar*(o1+o2)	
	ik = gk*(v-ek)		
} 

INITIAL{
	rate(celsius)
	SOLVE kin STEADYSTATE sparse
	
} 

KINETIC kin{ 
	acai =  (lcai)/diff 

	if (acai < cai)
		{acai = cai}
	rates(acai) 
	
	~c1<->c2 (dirc2_t_ca, invc1_t) 
	~c2<->c3 (dirc3_t_ca, invc2_t) 
	~c3<->c4 (dirc4_t_ca, invc3_t) 
	~c3<->o1 (diro1_t, invo1_t) 
	~c4<->o2 (diro2_t, invo2_t) 
	CONSERVE c1+c2+c3+c4+o2+o1=1 
} 

FUNCTION temper (Q10, celsius (degC)) {
	temper = Q10^((celsius -23(degC)) / 10(degC)) 
}

PROCEDURE rates(c(mM)){
	dirc2_t_ca = dirc2_t*c * fac
	dirc3_t_ca = dirc3_t*c * fac
	dirc4_t_ca = dirc4_t*c * fac
} 

PROCEDURE rate (celsius(degC)) {
	tcorr = temper (Q10,celsius)
	invc1_t = invc1*tcorr  
	invc2_t = invc2*tcorr
	invc3_t = invc3*tcorr 
	invo1_t = invo1*tcorr 
	invo2_t = invo2*tcorr 
	diro1_t = diro1*tcorr 
	diro2_t = diro2*tcorr 
	dirc2_t = dirc2*tcorr
	dirc3_t = dirc3*tcorr
	dirc4_t = dirc4*tcorr
}