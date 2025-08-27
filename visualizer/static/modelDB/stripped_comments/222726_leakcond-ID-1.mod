UNITS {
        (mV) = (millivolt)
        (mA) = (milliamp)
        (mM) = (milli/liter)
        (S) = (siemens)
        FARADAY = 96480 (coul)
        R       = 8.314 (volt-coul/degC)
}

NEURON {
        SUFFIX lc
        USEION na READ nai,nao WRITE ina
        USEION k READ ki,ko WRITE ik
        USEION ca READ cai,cao WRITE ica
	 
        RANGE ina,ik,ica,icl
        RANGE gna,gk,gca,gcl
        RANGE cnar,ckr,ccar,cclr
}

PARAMETER {
        v               (mV)    
        celsius         (degC)
        nai             (mM)
        nao             (mM)
        pna= 2.07e-7    (cm/s)
        gna=0.000005     (S/cm2)
        cnar=38.38      (mM)        
        ki              (mM)
        ko              (mM)
        pk= 3.45e-6     (cm/s)
        gk=0.00005      (S/cm2)
        ckr=15.77       (mM)        
        cai             (mM)
        cao             (mM)
        pca= 0          (cm/s)
        gca=0           (S/cm2)
        ccar=0.00184    (mM)        
	
	
	
	
	
}

ASSIGNED {
        ina     (mA/cm2)
        ik      (mA/cm2)
        ica     (mA/cm2)
	
        v0      (mV)
        vnar    (mV)
        vkr     (mV)
        vcar    (mV)
	vclr	(mV)
}

BREAKPOINT {
        v0 = (1000)*R*(celsius+273.16)/(FARADAY)
        vnar = v0*log(nao/nai)
        ina = gna * nao*nai*log(nao/nai)/(cnar*(nao-nai))* ghkfact(v,nai,nao,1)*(v-vnar)
        vkr = v0*log(ko/ki)
        ik = gk * ko*ki*log(ko/ki)/(ckr*(ko-ki))* ghkfact(v,ki,ko,1)*(v-vkr)
        vcar = v0*log(cao/cai)/2
        ica = gca * cao*cai*log(cao/cai)/(ccar*(cao-cai))*ghkfact(v,cai,cao,2)*(v-vcar)
	
	
}

FUNCTION expfun( v(mV), v0(mV)) (1) { 
        LOCAL e, w
        w=v/v0
        if (fabs(w)>1e-4) 
          { e = w / (exp(w)-1) }
        else 
        
          { e = 1-w/2+w^2/12 }
        expfun=e
}

FUNCTION ghkfact(v(mV), ci(mM), co(mM), z)  (1) { 
        LOCAL v0, vr
        v0 = (1000)*R*(celsius+273.16)/(z*FARADAY)
        vr = v0*log(co/ci)
        ghkfact = expfun(v,v0)/(expfun(v-vr,v0)*expfun(vr,v0))
}