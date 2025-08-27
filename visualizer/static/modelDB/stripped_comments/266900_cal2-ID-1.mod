UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)

	FARADAY = 96520 (coul)
	R = 8.3134 (joule/degC)
	KTOMV = .0853 (mV/degC)
}

PARAMETER { 
     
    curr
	v (mV)
	celsius 	(degC)
	gcalbar=.003 (mho/cm2)
	ki=.001 (mM)
	cai = 50.e-6 (mM)
	cao = 2 (mM)
	q10 = 5
	mmin=0.2
	tfa = 1
	a0m =0.1
	zetam = 2
	vhalfm = 4
	gmm=0.1	
	sh=0
	ggk
	time1=600
	time0=100
	alphash0=0
	alphash1=0.15
	sh2
	count=1
	vrun (mV)
	delta=0
	vinit=-76.2 
     vvrun=0	
	  timestep=1000
	 vrun2
	 v0
	dv0
	ddv
	flag=0
	FCa = 2
	PCa = 1
	BCa = 2
	CCa = 50
	stim_moltCa=1
	eca2=140
	
}


NEURON {
	SUFFIX caldb 
	POINTER stim_i
	USEION ca READ cai,cao WRITE ica
        RANGE flag, curr, gcalbar,cai, ica, gcal, ggk,sh,count,delta2,vrun2, stim_moltCa,eca2
        GLOBAL  minf,tau
}

STATE {
	m
}

ASSIGNED {
	ica (mA/cm2)
        gcal (mho/cm2)
        minf
        tau   (ms)
		stim_i
	
}

INITIAL {
	rate(v,sh2)
	m = minf
	vrun=0
	vvrun=vrun
}

BREAKPOINT {
	SOLVE state METHOD cnexp
	gcal = gcalbar*m*m*h2(cai)
	ggk=ghk(v,cai,cao)
	ica = gcal*ggk

}

FUNCTION h2(cai(mM)) {
	h2 = ki/(ki+cai)
}


FUNCTION ghk(v(mV), ci(mM), co(mM)) (mV) {
        LOCAL nu,f

        f = KTF(celsius)/2
        nu = v/f
        ghk=-f*(1. - (ci/co)*exp(nu))*efun(nu)
}

FUNCTION KTF(celsius (DegC)) (mV) {
        KTF = ((25./293.15)*(celsius + 273.15))
}


FUNCTION efun(z) {
	if (fabs(z) < 1e-4) {
		efun = 1 - z/2
	}else{
		efun = z/(exp(z) - 1)
	}
}

FUNCTION alp(v(mV)) (1/ms) {
	alp = 15.69*(-1.0*v+81.5)/(exp((-1.0*v+81.5)/10.0)-1.0)
}

FUNCTION bet(v(mV)) (1/ms) {
	bet = 0.29*exp(-v/10.86)
}

FUNCTION alpmt(v(mV)) {
  alpmt = exp(0.0378*zetam*(v-vhalfm)) 
}

FUNCTION betmt(v(mV)) {
  betmt = exp(0.0378*zetam*gmm*(v-vhalfm)) 
}

DERIVATIVE state {  
        rate(v,sh2)
        m' = (minf - m)/tau
}

BEFORE STEP { LOCAL i
       
	  if(stim_i==0 && flag==0){ 
		  vrun=0
		  vvrun=0
		  
	    }else{
		 flag=1
		             		  
		delta=v-vinit
		if (count<timestep+1){
		    vrun= (delta-vrun)*(FCa/(count+1))+vrun
			vrun2=vrun 
		 }else{

		vrun2= (delta)*(FCa/(timestep+1))+vrun2*pow((1-FCa/(timestep+1)),PCa)
			
			}
		 	     	   
	   vvrun=(BCa*vrun2/(1+vrun2/CCa))
	    
		count=count+1   
        }
		  						
		 sh2=sh+alphash1*vvrun
	
}
   

PROCEDURE rate(v (mV),sh2) { 
        LOCAL a, b, qt, i
        qt=q10^((celsius-25)/10)
			
        a = alp(v-sh2)
        b = 1/((a + bet(v-sh2)))
        minf = a*b
	tau = betmt(v-sh2)/(qt*a0m*(1+alpmt(v-sh2)))
	if (tau<mmin/qt) {tau=mmin/qt}
}