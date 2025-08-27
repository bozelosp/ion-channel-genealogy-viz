UNITS {
        (mA) = (milliamp)
        (mV) = (millivolt)
}

PARAMETER {
	celsius
        v (mV)
        gkabar=.008 (mho/cm2)
        vhalfn=-1   (mV)
        vhalfl=-56   (mV)
        a0l=0.05      (/ms)
        a0n=.1    (/ms)
        zetan=-1.8    (1)
        zetal=3    (1)
        gmn=0.39   (1)
        gml=1   (1)
        lmin=2  (mS)
        nmin=0.2  (mS)
        pw=-1    (1)
        tq=-40
        qq=5
        q10=5
        qtl=1
	ek
}


NEURON {
        SUFFIX kad
        USEION k READ ek WRITE ik
        RANGE gkabar,gka,vhalfn,vhalfl,i
        GLOBAL ninf,linf,taul,taun,lmin
}

STATE {
        n
        l
}

ASSIGNED {
        ik (mA/cm2)
        i (mA/cm2)
        ninf
        linf      
        taul
        taun
        gka
}

BREAKPOINT {
        SOLVE states METHOD cnexp
        gka = gkabar*n*l
        i = gka*(v-ek)
        ik = i

}

INITIAL {
	rates(v)
	n=ninf
	l=linf
}


FUNCTION alpn(v(mV)) {
    LOCAL zeta,Arg
    Arg=(v-tq)/qq
    
    if (Arg<-50) {zeta=zetan+pw}
    else if (Arg>50) {zeta=zetan}
    else {zeta=zetan+pw/(1+exp(Arg))}
    
    Arg=1.e-3*zeta*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))
    
    if (Arg<-50) {alpn=0}
    else if (Arg>50) {alpn=exp(50)}
    else {alpn=exp(Arg)}
}

FUNCTION betn(v(mV)) {
    LOCAL zeta,Arg
    Arg=(v-tq)/qq
    
    if (Arg<-50) {zeta=zetan+pw}
    else if (Arg>50) {zeta=zetan}
    else {zeta=zetan+pw/(1+exp(Arg))}
    
    Arg=1.e-3*zeta*gmn*(v-vhalfn)*9.648e4/(8.315*(273.16+celsius))
    
    if (Arg<-50) {betn=0}
    else if (Arg>50) {betn=exp(50)}
    else {betn=exp(Arg)}
}

FUNCTION alpl(v(mV)) {
    LOCAL Arg
    Arg=1.e-3*zetal*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))
    
    if (Arg<-50) {alpl=0}
    else if (Arg>50) {alpl=exp(50)}
    else {alpl=exp(Arg)}
}

FUNCTION betl(v(mV)) {
    LOCAL Arg
    Arg=1.e-3*zetal*gml*(v-vhalfl)*9.648e4/(8.315*(273.16+celsius))
    
    if (Arg<-50) {betl=0}
    else if (Arg>50) {betl=exp(50)}
    else {betl=exp(Arg)}
}

DERIVATIVE states {  
        rates(v)
        n' = (ninf - n)/taun
        l' = (linf - l)/taul
}

PROCEDURE rates(v (mV)) { 
        LOCAL a,qt
        qt=q10^((celsius-24)/10)
        a = alpn(v)
        ninf = 1/(1 + a)
        taun = betn(v)/(qt*a0n*(1+a))
        if (taun<nmin) {taun=nmin}
        a = alpl(v)
        linf = 1/(1+ a)
        taul = 0.26*(v+50)/qtl
        if (taul<lmin/qtl) {taul=lmin/qtl}
}