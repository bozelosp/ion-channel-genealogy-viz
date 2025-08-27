UNITS {
        (molar) = (1/liter)
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
         F = (faraday) (coulomb)
         R = (mole k) (mV-coulomb/degC)
        (mM) =  (millimolar)
}
 
NEURON {
        SUFFIX hhb
        
        USEION k WRITE ik
        RANGE  nspeed, sshift, shift,gnabar,gkhhbar,gkabar,ina,ikhh,ika,ik,ena,gk,gna,gka,nshift,bkcor, ashift, hshift, nslope, deptog, vhblock, taukv4,asshift,fptog, qfix, minf,hinf,qinf,pinf,ninf,hsinf,htau,mtau,ntau,hstau, qtau,ptau
        
}
 
INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}
 
PARAMETER {
        v   (mV)
        dt  (ms)
	nai (mM)
	celsius = 35.0 (degC)
        gnabar  = 550.0e-6 (S/cm2)
        gkhhbar = 665.0e-6 (S/cm2)
        gkabar  = 266.0e-6  (S/cm2) 
        ek  = -90.0  (mV)
        nao =  145  (mM)
        qv = 66.0 (mV)
        qs = 10.0  (1)
        shift=0 (mV)
        sshift=0 (mV)
        hshift = 0 (mV)
        ashift = 0.0 (mV)
        nshift = 0 (mV)
        nslope = 12
        bkcor = 0.0 (ms)
        deptog = 1 (1)
        vhblock = -54.0 (mV)
        taukv4 = 1.0
        asshift = 0
        fptog = 0
        qfix = 0
        nspeed = 2
 	
}
 
STATE {
         n <1e-4> p <1e-4> q <1e-4>
}
 
ASSIGNED {
        ina (mA/cm2)
        ik (mA/cm2)
        ika (mA/cm2)
        ikhh (mA/cm2)
        ena (mV)
	minf hinf ninf qinf pinf hsinf
	mtau htau hstau qtau ptau ntau
	
	gk 
	gka
}
 
BREAKPOINT {
        SOLVE states METHOD cnexp
	
    
    ena = 50
	
  	gk = n*n*n
	
	if(fptog){
	  gka = p*p*p*qfix
	} else {
	  gka = p*p*p*q 
	}
    
    ikhh = gkhhbar*gk*(v - ek)      
    ika = gkabar*gka*(v - ek)      
    ik = ika + ikhh

}
 
UNITSOFF
 




PROCEDURE rates(vf){
LOCAL ahf, bhf, ahs, bhs, am, bm
        minf = boltz(vf,-30.0-shift,9.7-sshift) 
        if(v > -54.0-hshift){
          hinf = boltz(vf,-54.0-hshift,-10.7) 
        } else {
          hinf = 0.5
        }
        hsinf = boltz(vf,vhblock,-1.57)
        
        
        

        ninf = boltz(vf,-20,nslope)

        
        
        
        pinf = boltz(vf,-35.0,7.0) 
        qinf = boltz(vf,-61.0+ashift,-4.5+asshift) 
        pinf = pow(pinf,1.0/3.0)
        
        if (-19.565-0.5*vf > 1e-3 || -19.565-0.5*vf < -1e-3){
          am = -(15.65 + 0.4*vf)/(exp(-19.5625-0.5*vf) - 1.0)
        } else {
          am = 0.8
        }
        bm = 3*exp(-7.5e-3*vf)
       
        mtau =0.01 + 1.0/(am+bm)
        
        
        
        ahf = 5.0e-4*exp(-6.3e-2*vf)
        bhf = 9.75*exp(0.134*vf)
        htau = 0.4 + 1.0/(ahf+bhf)
        ahs = 1.3e-2*exp(8.5e-4*vf)
        bhs = 13.165*exp(0.14556*vf)
        
        hstau = 20.0 +180.0*boltz(vf, 0,-1.0) 
										
										

        if (v > -40-nshift){
          ntau = 1 + 5.0*exp(-(log(1.0+0.05*(vf+40+nshift))/0.05*log(1.0+0.05*(vf+40+nshift))/0.05)/300) 
        } else {
          
          ntau = nspeed +(6-nspeed)*exp(-((vf+40+nshift))*((vf+40+nshift))/8.0) 
          
        }
        ntau =1.0*ntau
        ptau = 0.1*(1.029 + 4.83*boltz(vf,-56.7,-6.22))
        qtau = taukv4*250.0


}
INITIAL {
        rates(v)
        
        
        n = ninf
        p = pinf
        q = qinf
        
}

DERIVATIVE states {  
        rates(v)

        
        
        
        n' = 1.0*(ninf-n)/ntau
        p' = (pinf-p)/ptau
        q' = (qinf-q)/qtau
}
 
 
 
FUNCTION boltz(x,y,z) {
                boltz = 1/(1 + exp(-(x - y)/z))
}
 
UNITSON