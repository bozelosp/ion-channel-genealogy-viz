NEURON {
	POINT_PROCESS GrCAMPAplusNMDA
	NONSPECIFIC_CURRENT i
	RANGE Q10_diff,Q10_channel, fix_celsius
	RANGE g_nmda, g_ampa, i_ampa, i_nmda
	RANGE Cdur, e_nmda, e_ampa
	RANGE r1FIX, r2, r3,r4,r5, gmax_ampa, r1,r6,r6FIX,kB
	RANGE Rb, Ru, Rr, Ro, Rc, rb, gmax_nmda, RdRate
	RANGE tau_1, tau_rec, tau_facil, U
	RANGE PRE,T,Tmax
	RANGE Trelease, diffuse
  RANGE Trscale
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(umho) = (micromho)
	(mM) = (milli/liter)
	(pS) = (picosiemens)
	(nS) = (nanosiemens)
	(uS) = (microsiemens)
	(um) = (micrometer)
	PI	= (pi)		(1)
  }

  PARAMETER {
	Q10_diff	= 1.5
	Q10_channel	= 2.4
	
	gmax_ampa		= 1200e-6   (uS)
	gmax_nmda		= 18800e-6  	(uS)	

	
	r1FIX		= 5.4		(/ms/mM)
	r2			= 0.82	(/ms)
	r3			= 0		(/ms)
	r4			= 0		(/ms)
	r5			= 0.013	(/ms)
	r6FIX		= 1.12	(/ms/mM)
	e_ampa		= 0	(mV)
	kB			= 0.44	(mM)

	
	Rb	=  5		(/ms/mM)  	
	Ru	=  0.1		(/ms)		
	RdRate	=  12e-4  	(/ms)		
	Rr	=  9e-3		(/ms)		
	Ro	=  3e-2 	(/ms)		
	Rc	=  0.966	(/ms)		
	e_nmda	= -3.7  (mV)	
	mgo = 1 (mM)
	kmg = 1.77 (mM)
	deltav = 22.4 (mV)


	
	Cdur		= 0.3	(ms)
	U 			= 0.42 (1) 	< 0, 1 >
	tau_rec 	= 8 (ms) 	< 1e-9, 1e9 >
	tau_facil 	= 5 (ms) 	< 0, 1e9 >
	tau_1 		= 1 (ms) 	< 1e-9, 1e9 >

	u0 			= 0 (1) 	< 0, 1 >	
	Tmax		= 1  (mM)

  Rd	= 1.03 (um)
	Diff	= 0.223 (um2/ms)
  Trscale = 0.83 (mM)

	fix_celsius = 37 (degC)

}


ASSIGNED {
	v		(mV)		
	i 		(nA)		

	g_ampa 		(pS)		
	i_ampa  (nA)

	r1		(/ms)
	r6		(/ms)
	T		(mM)
	Trelease	(mM)
	diffuse (mM)

	g_nmda    (pS)
	i_nmda  (nA)
	rb		(/ms)    
	MgBlock

	tsyn		(ms)
	tr1

	gbar_Q10 (mho/cm2)
	Q10 (1)
  r_fast (/ms)
  r_slow (/ms)

  gmx_ampa
	gmx_nmda

  df11 (/ms)
  df12 (/ms)
  df21 (/ms)
  df22 (/ms)
  df23 (/ms)
  df31 (/ms)
  df32 (/ms)
  df33 (/ms)
  df41 (/ms)
}

STATE {
	O_ampa
	D_ampa

	
	C1_nmda		
	C2_nmda		
	 D_nmda		
	 O_nmda		

	PRE
	Tr
	I1
	I2
	I3
}

INITIAL {
	O_ampa=0
	D_ampa=0
	T=0 (mM)

	C1_nmda = 0
	C2_nmda = 0
	D_nmda  = 0
	O_nmda  = 0

	PRE = 0
	Tr = 0
	I1 = 0
	I2 = 0
  I3 = 0

	Trelease=0 (mM)
	diffuse=0 (mM)


	gbar_Q10 = Q10_diff^((fix_celsius-30)/10)
	Q10 = Q10_channel^((fix_celsius-30)/10)

	if(tau_1>=tau_rec){
		printf("Warning
		tau_rec=tau_1+1e-5
	}

    r_fast = 4/(Rd*Rd/(4*Diff))
    r_slow = 0.01/(Rd*Rd/(4*Diff))

    df11 = -6.2*r_slow -20*r_slow
    df12 = 20*r_slow
    df21 = 9.09*r_slow
    df22 = -4.9*r_slow
    df23 = df22-df21
    df31 = 1.71*r_slow
    df32 = -0.55*r_slow
    df33 = df32-df31
    df41 = 0.333*r_slow

    gmx_ampa = gmax_ampa * gbar_Q10
		gmx_nmda = gmax_nmda * gbar_Q10

		g_ampa = 0
		g_nmda = 0

}

BREAKPOINT {
		LOCAL tact
		diffuse = Tr * Trscale
    Trelease = T+diffuse
    tact = Trelease*Trelease/(Trelease + kB)/(Trelease + kB)
    r1 = r1FIX * tact
    r6 = r6FIX * tact

		rb = Rb * diffuse

		rates(v)

    SOLVE transition METHOD cnexp

    g_ampa = gmx_ampa * O_ampa
		g_nmda = gmx_nmda * O_nmda * MgBlock
		i_ampa = g_ampa * (v - e_ampa)
		i_nmda = g_nmda * (v - e_nmda)
    i =  i_ampa + i_nmda
}

DERIVATIVE transition {
    PRE' = -r_fast*PRE
    Tr' = df11*Tr + PRE + df12*I1
    I1' = df21*Tr + df23*I1 - df22*I2
    I2' = df31*I1 + df33*I2 - df32*I3
    I3' = df41*(I2-I3)

    O_ampa' = (r1*(1-O_ampa-D_ampa) - r2*O_ampa)*Q10
    D_ampa' = (-r5*D_ampa + r6*(1-O_ampa-D_ampa))*Q10

		C1_nmda' = (rb*(1-C1_nmda-C2_nmda-D_nmda-O_nmda) + Ru*C2_nmda - (rb+Ru)*C1_nmda)*Q10
		C2_nmda' = (rb*C1_nmda + Rr*D_nmda - (Ru+RdRate)*C2_nmda)*Q10
		D_nmda'  = (RdRate*C2_nmda - Rr*D_nmda)*Q10
		O_nmda'  = (Ro*C2_nmda - Rc*O_nmda)*Q10

}


PROCEDURE rates(v(mV)) {
	
	TABLE MgBlock DEPEND mgo,kmg,deltav FROM -120 TO 30 WITH 150
	MgBlock = 1 / ( 1 + (mgo/kmg)*exp ( - v/deltav ) )
}


NET_RECEIVE(weight, on, x, y, z, u, tsyn (ms)) {LOCAL efr, ef1, dtsyn







INITIAL {
    x = 0
	y = 0
	z = 0
	u = u0
	tsyn = t
  tr1 = 1/((tau_1/tau_rec)-1)
}
   if (flag == 0) {
		if (!on) {
			on = 1
            dtsyn = t - tsyn
            efr = exp(-dtsyn/tau_rec)
            ef1 = exp(-dtsyn/tau_1)

			z = z*efr
			z = z + y*(ef1 - efr)*tr1
			y = y*ef1
			x = 1-y-z

			if (tau_facil > 0) {
				u = u*exp(-dtsyn/tau_facil)
				u = u + U * ( 1 - u )
			} else { u = U }
			y = y + x * u

			T=Tmax*y
			PRE = PRE + y

			tsyn = t
			net_send(Cdur, 1)
		} else {
			net_move(t+Cdur)
		}
  }
  if (flag == 1) {
		T = 0
		on = 0
  }
}