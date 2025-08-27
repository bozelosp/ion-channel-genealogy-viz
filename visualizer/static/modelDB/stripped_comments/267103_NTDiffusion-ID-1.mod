DEFINE NDates 10


NEURON {
	POINT_PROCESS NTDiffusion
	RANGE Radius, CleftWidth, Diffusivity, BasalNTConcentration, NTi, k, Nused
	RANGE NTConcentration
	RANGE comp
	RANGE NTRatio
	RANGE tDiff
}

UNITS {
	(molar) = (1/liter)			
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
	Radius	= 1	(um)			
	CleftWidth	= 0.02	(um)		
	Diffusivity	= 0.33 (um.um/ms)	
	BasalNTConcentration = 0.0 (mM)	
	k = 1.0 					
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
			
			NTRatio = 0
		}
	
		if(Radius == 0){
		comp = NTRatio * 1e18 * NTi / (4 * PI * CleftWidth * Diffusivity) 
		} else {
		comp = NTRatio * 1e18 * NTi / (4 * PI * CleftWidth * Diffusivity) * exp(- NTRatio * (Radius*Radius)/(4*Diffusivity)) 
		}
		
		
		if(comp >= 0.001*NTConcentration){
			NTConcentration = NTConcentration + comp
		} else if(t==tr[i]){
			NTConcentration = BasalNTConcentration
		} else {
			Nused = i
			
		}
	}
	
	
	
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

	if(flag == 0){ 
		if(Nused < NDates){
			Nused = Nused +1
		} else {
			printf("Size of the NTDiffusion Glu array might be too small")
		}
		
		NTi = k*3000/Nav 
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
		
	}
}