NEURON	{
	SUFFIX TC_ih_Bud97
	NONSPECIFIC_CURRENT ih
	RANGE gh_max, g_h, i_rec
        RANGE shift, vhalf, k, a, b, tvhalf1, tvhalf2, mInf, m
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gh_max = 2.2e-5 (S/cm2) 
	e_h =  -43.0 (mV)
        celsius (degC)
	q10 = 4 
        shift=11
        vhalf=-86.4
        k=11.2
        a=0.086
        b=0.0701
        tvhalf1=-169.6511627906977
        tvhalf2=26.67617689015692		     
}

ASSIGNED	{
	v	(mV)
	ih	(mA/cm2)
	g_h	(S/cm2)
	mInf
	mTau
	tcorr		
	i_rec
}

STATE	{ 
	m
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	g_h = gh_max*m
	ih = g_h*(v-e_h)
	i_rec = ih
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
}

INITIAL{
  	tcorr = q10^((celsius-34)/10)  
	rates()
	m = mInf
}

UNITSOFF
PROCEDURE rates(){
        mInf = 1/(1+exp((v-vhalf+shift)/k)) 
        mTau = (1/(exp(-(-tvhalf1 + v+shift)*a) + exp((-tvhalf2 + v+shift)*b )))/tcorr 
}
UNITSON