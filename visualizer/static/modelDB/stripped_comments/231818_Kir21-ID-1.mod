NEURON {
	SUFFIX Kir21
	USEION k READ ek WRITE ik
    RANGE  ik, gk, gkbar
	GLOBAL mg_i, As, shiftmg, cas,fac, gsub, b, spm_i, vshiftbs
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
        (S)  = (siemens)
	
	(molar) = (1/liter)
	(mM) = (millimolar)
	(uM) = (micromolar)
	
}

PARAMETER {
	v 		(mV)
	gkbar  = 0.00015                (S/cm2) 
	mg_i = 4 (mM)  
	spm_i = 1 (uM) 
	As = 1
	vshiftbs = 0 (mV)
	b= 0.1  
	
	
	fac = 0.0005  
	gsub = 0.1  
	shiftmg = 1 
	cas = 1  
}

STATE {
        O BS B1 B2 B3 BB
}

ASSIGNED {
        
        
        ik                             (mA/cm2)
        gk                            (S/cm2)
        ek                            (mV)
		alpha1   					(/ms)
		beta1						(/ms)
		alphas   					(/ms)
		betas   					(/ms)
		alphas2 					(/ms)
		betas2						(/ms)
}

INITIAL {
	rate(v)
	SOLVE kin STEADYSTATE sparse
}

BREAKPOINT {
	SOLVE kin METHOD sparse
	gk = (gkbar ) * (O + 1/3 * B2 + 2/3 * B1) + (gkbar * gsub ) * BS 
    ik = gk * ( v - ek )
}


KINETIC kin {
LOCAL alpha2, alpha3, beta2, beta3
rate(v)

alpha2 = 2*alpha1
beta2 = 2 * beta1
alpha3 = 3*alpha1
beta3 = 3*beta1


~ BS <-> O (alphas,betas)
~ B1 <-> O (alpha1,beta3)
~ B2 <-> B1 (alpha2,beta2)
~ B3 <-> B2 (alpha3,beta1)
~ BB <-> BS (alphas2,betas2)

CONSERVE O + BS + BB + B1 + B2 + B3 = 1

}


PROCEDURE rate(v (mV)) { 
	LOCAL a,d
	
	
	alpha1 = 12 * exp(-0.025 * (v - (shiftmg * (ek))))			
	beta1 = mg_i/8 * 28 * exp(0.025 * (v - (shiftmg * (ek))) ) 
	
	
	alphas = As * 0.17 * exp(cas*-0.07 * (v - (ek) +8/8 (mV/mM) * mg_i)) 
	
	betas =  As * spm_i * 0.28 * exp(0.15 * (v - (ek) +8/8  (mV/mM) * mg_i)) 
	
	
	

	a = - 1/9.1 + b
	
	

	
	alphas2 = fac* 40 * exp(a*(v-(ek+vshiftbs)))  
	betas2 = spm_i * fac * exp(b*(v-(ek+vshiftbs))) 

}