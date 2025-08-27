NEURON{
	SUFFIX Golgi_SK2
	USEION ca READ cai
	USEION k READ ek WRITE ik
	RANGE Q10_diff, Q10_gbar,Q10_channel,gbar_Q10, diff_Q10, fix_celsius
	RANGE gbar, g, ik, tcorr, diff, ic, o ,c
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) = (millimolar)
}

PARAMETER {
    fix_celsius = 37 (degC)
	cai (mM)
	gbar = 0.045 (mho/cm2)
	Q10_gbar	= 1.5
	Q10_diff	= 1.3
	Q10_channel	= 3
	diff = .333 (1)


	invc1 = 80e-3  ( /ms)
	invc2 = 80e-3  ( /ms)
	invc3 = 200e-3 ( /ms)

	invo1 = 1      ( /ms)
	invo2 = 100e-3 ( /ms)
	diro1 = 160e-3 ( /ms)
	diro2 = 1.2    ( /ms)


	dirc2 = 200 ( /ms-mM )
	dirc3 = 160 ( /ms-mM )
	dirc4 = 80  ( /ms-mM )
}

ASSIGNED{
	v	(mV)
	ek	(mV)
	g	(mho/cm2)
	ik	(mA/cm2)
	invc1_t  ( /ms)
	invc2_t  ( /ms)
	invc3_t  ( /ms)
	invo1_t  ( /ms)
	invo2_t  ( /ms)
	diro1_t  ( /ms)
	diro2_t  ( /ms)
	dirc2_t  ( /ms)
	dirc3_t  ( /ms)
	dirc4_t  ( /ms)
	tcorr	 (1)

	dirc2_t_ca  ( /ms-mM)
	dirc3_t_ca  ( /ms-mM)
	dirc4_t_ca  ( /ms-mM)
	gbar_Q10 (mho/cm2)
	diff_Q10 (mho/cm2)
	ic
	o
	c
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
	g = gbar_Q10*(o1+o2)	
	ik = g*(v-ek)		
	ic = ik
	o = o1+o2
	c = c1+c2+c3+c4
}

INITIAL{
	gbar_Q10 = gbar * (Q10_gbar^((fix_celsius-23)/10))
	diff_Q10 = diff * (Q10_diff^((fix_celsius-23)/10))
	rate(fix_celsius)
	SOLVE kin STEADYSTATE sparse
}

KINETIC kin{
	rates(cai*diff_Q10)
	~c1<->c2 (dirc2_t_ca, invc1_t)
	~c2<->c3 (dirc3_t_ca, invc2_t)
	~c3<->c4 (dirc4_t_ca, invc3_t)
	~c3<->o1 (diro1_t, invo1_t)
	~c4<->o2 (diro2_t, invo2_t)
	CONSERVE c1+c2+c3+c4+o2+o1=1
}

FUNCTION temper (Q10, fix_celsius (degC)) {
	temper = Q10^((fix_celsius -23(degC)) / 10(degC))
}

PROCEDURE rates(cai(mM)){
	dirc2_t_ca = dirc2_t*cai
	dirc3_t_ca = dirc3_t*cai
	dirc4_t_ca = dirc4_t*cai
}

PROCEDURE rate (fix_celsius(degC)) {
	tcorr = temper (Q10_channel,fix_celsius)
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