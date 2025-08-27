NEURON {
	POINT_PROCESS ExpSynSTDP_triplet
	RANGE tau, e, i, factor, gw
	RANGE R1, R2, O1, O2
	RANGE A2mais, A2menos, A3mais, A3menos
	RANGE Taumais, Taumenos, Tauy, Taux
	RANGE t_last_pre, t_last_post 
	RANGE flag_pre, flag_post 
	RANGE flag1_pre, flag1_post 
	RANGE peso
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

PARAMETER {
	e = 0	(mV)
	tau= 10 (ms)
	factor = 1
	gw=0.01
	latency=0
	A2mais=4.6e-3 
	A2menos=3e-3
	A3mais=9.1e-3
	A3menos=7.5e-9
	Taumais=16.8
	Taumenos=33.7
	Tauy=47
	Taux=575

	}

ASSIGNED {
	v (mV)
	i (nA)
	tpost (ms)
	t_last_pre (ms)
	t_last_post (ms)
	flag_pre
	flag_post
	flag1_pre
	flag1_post

}

STATE {
	g (uS)
	R1
	R2
	O1
	O2
	peso
}

INITIAL {
	g=0
	peso=0
	tpost = -1e9
	net_send(0, 1)
	O1=0
	O2=0
	R1=0
	R2=0	
	t_last_pre=0
	t_last_post=0
	flag_pre=0
	flag_post=0
	flag1_pre=0
	flag1_post=0
}

BREAKPOINT {
	SOLVE state METHOD cnexp

	if (g<0) {g=0} 


	
if (flag1_pre==1)  {
	peso= peso -factor*O1*(A2menos+A3menos*R2)

}

if (flag1_post==1)  {

	peso= peso + factor*R1*(A2mais*(1-peso/(1e3*gw))+A3mais*(1-peso/(1e3*gw))*O2) 
}

	i = g*(v - e)
}

DERIVATIVE state {




if (flag1_pre==0) {
	R1'= -R1/Taumais/1
	R2'= -R2/Taux/1
} 
if (flag1_pre==1)  {
	R1'= (-R1 + 1)/Taumais/1
	R2'= (-R2 + 1)/Taux/1
	flag1_pre=0 

}

if (flag1_post==0) {
	O1'= -O1/Taumenos/1
	O2'= -O2/Tauy/1
}
if(flag1_post==1)  {
	O1'= (-O1 + 1)/Taumenos/1
	O2'= (-O2 + 1)/Tauy/1
	flag1_post=0 

	}	

	g'= -g/tau


}

NET_RECEIVE(w (uS), A, tpre (ms)) {
	INITIAL { A = 0  tpre = 0 }


	if (flag == 0) { 

g = g + w + peso

if (g<0) {g=0} 

		tpre = t
		A = A 
		flag_pre=1 


	}else if (flag == 2) { 
		tpost = t

		FOR_NETCONS(w1, A1, tp) { 
			A1 = A1
		flag_post=1 
		flag_pre=0 
		}

if (flag_post ==1) { 

	flag1_post=1 
	t_last_post=t
	}


	} else { 
		WATCH (v > -20) 2
	}

if (flag_pre ==1) { 
	t_last_pre=t
	flag1_pre=1 

	}

}