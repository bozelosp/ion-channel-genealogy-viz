NEURON	{
	SUFFIX Kv11
	USEION k READ ek WRITE ik
	RANGE gkbar, gk, ik 
	GLOBAL md,mk,mtA,mtd,mtk
}

UNITS	{
	(S) = (siemens)
	(mV) = (millivolt)
	(mA) = (milliamp)
}

PARAMETER	{
	gkbar = 0.00001 (S/cm2) 
	md = -30.5 (mV)
	mk = -11.3943 (mV)
	mtA = 30 (ms)
	mtd = -76.56 (mV)
	mtk = 26.1479 (mV)
}

ASSIGNED	{
	v	(mV)
	ek	(mV)
	ik	(mA/cm2)
	gk	(S/cm2)
	mInf
	mTau
	hInf
	hTau
}

STATE	{ 
	m
	h
}

BREAKPOINT	{
	SOLVE states METHOD cnexp
	gk = gkbar*m*h
	ik = gk*(v-ek)
}

DERIVATIVE states	{
	rates()
	m' = (mInf-m)/mTau
	h' = (hInf-h)/hTau
}

INITIAL{
	rates()
	m = mInf
	h = hInf
}

PROCEDURE rates(){
	UNITSOFF 
		mInf = 1.0000/(1+ exp((v - md)/mk)) 
		mTau = mtA/(1+ exp((v - mtd)/mtk)) 
		hInf = 1.0000/(1+ exp((v - -30.0000)/27.3943)) 
		hTau = 15000.0000/(1+ exp((v - -160.5600)/-100.0000))
	UNITSON
}