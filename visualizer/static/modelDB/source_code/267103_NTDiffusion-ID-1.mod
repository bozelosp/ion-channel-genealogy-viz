COMMENT
ALE - 08.06.2011
NTDiffusion model using analytic solution [A1.2] of savtchenko - The optimal height of the synaptic cleft - 2007
ENDCOMMENT

DEFINE NDates 10
:DEFINE NDates 20

NEURON {
	POINT_PROCESS NTDiffusion
	RANGE Radius, CleftWidth, Diffusivity, BasalNTConcentration, NTi, k, Nused
	RANGE NTConcentration
	RANGE comp
	RANGE NTRatio
	RANGE tDiff
}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
}

CONSTANT{
	PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062862089986280348253421170679
	Nav = 6.02214179e23
}

PARAMETER {
	Radius	= 1	(um)			: Radius, position of the receptor (DEFAULT 1) (50 nanometers is nominal)
	CleftWidth	= 0.02	(um)		: CleftWidth
	Diffusivity	= 0.33 (um.um/ms)	: Diffusivity coefficient of the NT in the cleft (DEFAULT 0.4)
	BasalNTConcentration = 0.0 (mM)	: resting NT concentration
	k = 1.0 					: variable d'ajustement
	comp = 0
	tRound = 0
	tDiff = 0
	NTRatio = 0
}

ASSIGNED {
	NTConcentration
	NTi
	tr[NDates]
	Nused
}

STATE {
	NT1
	NT2
	NT3
	NT4
	NT5
	NT6
	NT7
	NT8
	NT9
}

INITIAL {
	NTConcentration = BasalNTConcentration
	NTi = 0
	FROM i = 0 TO NDates-1{
		tr[i] = -1e12
	}
	Nused = 0

}

: LOCAL comp

BREAKPOINT {
	SOLVE state METHOD cnexp
	NTConcentration = BasalNTConcentration
	
	FROM i = 0 TO Nused-1{
		if(i == 0) {
			NTRatio = NT1
			tDiff = t-tr[i]
		} else if(i == 1) {
			NTRatio = NT2
		} else if(i == 2) {
			NTRatio = NT3
		} else if(i == 3) {
			NTRatio = NT4
		} else if(i == 4) {
			NTRatio = NT5
		} else if(i == 5) {
			NTRatio = NT6
		} else if(i == 6) {
			NTRatio = NT7
		} else if(i == 7) {
			NTRatio = NT8
		} else if(i == 8) {
			NTRatio = NT9
		} else {
			:printf("Warning: More NT events than expected")
			NTRatio = 0
		}
	
		if(Radius == 0){
		comp = NTRatio * 1e18 * NTi / (4 * PI * CleftWidth * Diffusivity) : en mM
		} else {
		comp = NTRatio * 1e18 * NTi / (4 * PI * CleftWidth * Diffusivity) * exp(- NTRatio * (Radius*Radius)/(4*Diffusivity)) : en mM
		}
		:printf("comp : %4.5f\n", comp)
		:printf("t = %f :\t %1.0f / %1.0f\t : comp = %f \n", t, i+1, Nused, comp)
		if(comp >= 0.001*NTConcentration){
			NTConcentration = NTConcentration + comp
		} else if(t==tr[i]){
			NTConcentration = BasalNTConcentration
		} else {
			Nused = i
			:printf("component used until %1.0f, %f / %f \n", Nused, comp, NTConcentration)
		}
	}
	
	: NTConcentration = NTConcentration/.06 : testing modifying the amplitude of glutamate
	
}

DERIVATIVE state {

	NT1' = - (NT1 * NT1)
	NT2' = - (NT2 * NT2)
	NT3' = - (NT3 * NT3)
	NT4' = - (NT4 * NT4)
	NT5' = - (NT5 * NT5)
	NT6' = - (NT6 * NT6)
	NT7' = - (NT7 * NT7)
	NT8' = - (NT8 * NT8)
	NT9' = - (NT9 * NT9)

}

NET_RECEIVE(w){

	if(flag == 0){ : external event coming from netcon object
		if(Nused < NDates){
			Nused = Nused +1
		} else {
			printf("Size of the NTDiffusion Glu array might be too small")
		}
		
		NTi = k*3000/Nav : en moles :  (Nav*PI*1e-9*ri*ri*CleftWidth)
		FROM i = 1 TO Nused{
			tr[Nused+1-i] = tr[Nused-i]
			
		}

		tr[0] = t
		NT9 = NT8
		NT8 = NT7
		NT7 = NT6
		NT6 = NT5
		NT5 = NT4
		NT4 = NT3
		NT3 = NT2
		NT2 = NT1
		NT1 = 1
		:printf("NTi : %4.5f\n", NTi)
	}
}

