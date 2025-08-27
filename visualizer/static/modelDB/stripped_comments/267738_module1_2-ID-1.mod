UNITS { }

NEURON {
	SUFFIX CaSP
	
	
	RANGE K1, K2, K3, K4, K5i, K6i, K  
	RANGE Rmax, Umax, t1, t2
	RANGE phi1, phi2, phi3, phi4

	
	RANGE C1i, C1n1, C1n2, C1n3, C1n4 
	RANGE C2i, C2n1, C2n2, C2n3, C2n4
	RANGE C3, C4, C5, alpha

	USEION mg WRITE mgi VALENCE 2
	USEION cl READ cli
}

PARAMETER {
	
	K1 = 3000		
	K2 = 3			
	K3 = 400		
	K4 = 1			
	K5i = 4e5		
	K6i = 150		
	K = 850			
	SF_AM = 5
	Rmax = 10		
	Umax = 2000		
	t1 = 1			
	t2 = 13			
	CS0 = 0.03     	
	B0 = 0.00043	
	T0 = 0.00007 	
	phi1 = 0.004
	phi2 = 0.98
	phi3 = 0.0002
	phi4 = 0.999

	
	C1i = 0.154
	C1n1 = 0.01
	C1n2 = 0.15
	C1n3 = 0.01
	C1n4 = 85	
	C2i = 0.11
	C2n1 = -0.0315
	C2n2 = 0.27
	C2n3 = 0.015
	C2n4 = 70	
	C3 = 54.717
	C4 = -18.847
	C5 = 3.905
	alpha = 1.65

	
	vth = -55	
	spK_index=0
}

STATE {
	CaSR			
	CaSRCS			
	Ca				
	CaP				
	CaT				
	AM				
	C1
	C2
	mgi				
}

ASSIGNED {
	v 	(mV)
	R
	t_shift 		
	R_On 			
	SpiKe_On 		
	K5
	K6
	AMinf
	AMtau
	C1inf
	C1tau
	C2inf
	C2tau
	cli			
	spK[100] 	
	xm[2] 		
	vm
	phi0
}

BREAKPOINT { LOCAL i, temp_R
	
	SPK_DETECT (v, t) 
	CaR (CaSR, t)

	SOLVE state METHOD cnexp
	
	mgi = AM^alpha

	xm[0]=xm[1]
	xm[1]=cli
	vm = (xm[1]-xm[0])/(dt*10^-3)
}

DERIVATIVE state {
	rate (cli, CaT, AM, t)
	
	CaSR' = -K1*CS0*CaSR + (K1*CaSR+K2)*CaSRCS - R + U(Ca)
	CaSRCS' = K1*CS0*CaSR - (K1*CaSR+K2)*CaSRCS
	
	Ca' = - K5*T0*Ca + (K5*Ca+K6)*CaT - K3*B0*Ca + (K3*Ca+K4)*CaP + R - U(Ca)
	CaP' = K3*B0*Ca - (K3*Ca+K4)*CaP
	CaT' = K5*T0*Ca - (K5*Ca+K6)*CaT
	
	AM' = (AMinf - AM)/AMtau
	C1' = (C1inf - C1)/C1tau
	C2' = (C2inf - C2)/C2tau
	mgi' = 0
}

PROCEDURE SPK_DETECT (v (mv), t (ms)) {
	if (SpiKe_On == 0 && v > vth) {
	SpiKe_On = 1
	spK[spK_index] = t
	spK_index = spK_index + 1
	R_On = 1
	} else if (v < vth) {
	SpiKe_On = 0
	}
}

FUNCTION U (x) {
	if (x >= 0) {U = Umax*(x^2*K^2/(1+x*K+x^2*K^2))^2}
	else {U = 0}
}

FUNCTION phi (x) {
	if (x <= 5) {
	phi = phi1*x + phi2
	phi0 = phi}
	else {
	phi = phi3*x + phi4
	phi0 = phi}
}

PROCEDURE CaR (CaSR (M), t (ms)) { LOCAL i, temp_R  
	if (R_On == 1) {
		FROM i=0 TO spK_index-1 {
			temp_R = temp_R + CaSR*Rmax*(1-exp(-(t-spK[i])/t1))*exp(-(t-spK[i])/t2)
		}
		R = temp_R
		temp_R = 0
	}
	else {R = 0}
}

PROCEDURE rate (cli (M), CaT (M), AM (M), t(ms)) {
	K5 = phi(cli)*K5i
	K6 = K6i/(1 + SF_AM*AM)
	AMinf = 0.5*(1+tanh(((CaT/T0)-C1)/C2))
	AMtau = C3/(cosh(((CaT/T0)-C4)/(2*C5)))
	C1inf = C1n1*(1+tanh(((CaT/T0)-C1n2)/C1n3)) + C1i
	C1tau = C1n4
	C2inf = C2n1*(1+tanh(((CaT/T0)-C2n2)/C2n3)) + C2i
	C2tau = C2n4
}

INITIAL {LOCAL i
	CaSR = 0.0025  	
	CaSRCS = 0		
	Ca = 1e-10		
	CaP = 0			
	CaT = 0			
	AM = 0			
	C1=C1i
	C2=C2i
	mgi = 0
		
	FROM i = 0 TO 99 {
	spK[i] = 0
	}
	spK_index = 0
	R_On = 0
	phi0 = 0
}