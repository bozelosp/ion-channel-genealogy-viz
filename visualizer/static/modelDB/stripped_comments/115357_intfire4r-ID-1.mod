NEURON {
	ARTIFICIAL_CELL IntFire4r
	RANGE taue, taui1, taui2, taum, e, i1, i2, m, ae, ai2, ai1
	RANGE nself, nexcite, ninhibit
	GLOBAL eps, taueps
}

PARAMETER {
	taue = 5 (ms)      <0,1e9>
	taui1 = 10 (ms)    <0,1e9>
	taui2 = 20 (ms)    <0,1e9>
	taum = 50 (ms)     <0,1e9>
	ib = 0
	eps = 1e-6         <1e-9,0.01>
	taueps = 0.01      <1e-9,1e9>
}

ASSIGNED {
	e i1 i2 m
	enew i1new i2new mnew
	t0 (ms)
	nself nexcite ninhibit

	ke (1/ms) ki1 (1/ms) ki2 (1/ms) km (1/ms)
	ae (1/ms) ai1 (1/ms) ai2 (1/ms)
	be bi1 bi2
	a b
	tau_swap (ms)
	flag (1)
}

PROCEDURE newstates(d(ms)) {
	LOCAL ee, ei1, ei2, em
	
	ee = exp(-ke*d)
	ei1 = exp(-ki1*d)
	ei2= exp(-ki2*d)
	em = exp(-km*d)
	enew = e*ee
	i1new = i1*ei1
	i2new = i2*ei2 + bi1*i1*(ei2 - ei1)
	mnew = m*em
		+ be*e*(em - ee)
		+ (bi2*i2 + a*i1)*(em - ei2)
		- b*i1*(em - ei1)
}

FUNCTION M() {
	newstates(t - t0)
	M = mnew
}

FUNCTION E() {
	newstates(t - t0)
	E = ae*enew
}

FUNCTION I() {
	newstates(t - t0)
	I = ai2*i2new
}

PROCEDURE update() {
	e = enew
	i1 = i1new
	i2 = i2new
	m = mnew
	t0 = t
}

PROCEDURE fixprecondition() {
    
	
	
    
    
    
    
    
    
    
    
    
    
    if(taui2 < 4*taueps){
        taui2=4*taueps
    }
    if(taui1 < 3*taueps){
        taui1=3*taueps
    }
    
    if(taue < 2*taueps){
        taue=2*taueps
    }
    
    
    if (taue > taui2) { 
        tau_swap=taue
        taue = taui2 - taueps 
        printf("Warning
    } else if (taui2-taue < taueps){
        taue = taui2 - taueps 
    }

    
	if (taui1 > taui2) {
        tau_swap=taui2
        taui2=taui1
        taui1=tau_swap
        printf("Warning
    }   
    
    
    
    if (taui2-taui1 < taueps){          
        taui1 = taui2 - taueps 
    }

    if (taum <= taui2) {
        if (taui2 -taum < taueps){      
            taum=taui2-taueps
        }
    
        if (fabs(taui1 -taum) < taueps){
            taum=taui1-taueps
        }
    
        if (fabs(taui1 -taum) < taueps){
            if(taui1 -taum < 0){
                taum=taui1-taueps
            } else{
                taui1=taum-taueps
            }
        }
    
        if (fabs(taue -taum) < taueps){
            if(taue -taum < 0){
                taum=taue-taueps
            } else{
                taue=taum-taueps
            }
        }
    
        if (fabs(taui1 -taum) < taueps){
            taum=taui1-taueps
        }
    
    } else if (taum-taui2 < taueps){    
        taum=taui2+taueps
    }
}

PROCEDURE factors() {
	LOCAL tp
    
	ke=1/taue  ki1=1/taui1  ki2=1/taui2  km=1/taum

	
	tp = log(km/ke)/(km - ke)
	be = 1/(exp(-km*tp) - exp(-ke*tp))
	ae = be*(ke - km)
        
	
    
	tp = log(ki2/ki1)/(ki2 - ki1)
	bi1 = 1/(exp(-ki2*tp) - exp(-ki1*tp))
	ai1 = bi1*(ki1 - ki2)
    
	
    
	
	e = 0
	i1 = 1
	i2 = 0
	m = 0
	bi2 = 1
	ai2 = bi2*(ki2 - km)
	a = bi2*bi1
	b = a*(ki2 - km)/(ki1 - km)
	
    
	tp = search()
	
    
	newstates(tp)
	bi2 = 1/mnew
	ai2 = bi2*(ki2 - km)
	a = bi2*bi1
	b = a*(ki2 - km)/(ki1 - km)
	
	
	newstates(tp)
    
	i1 = 0
}

FUNCTION deriv(d(ms)) (/ms2) { 
    deriv= ( - km *exp( - d*km ) + ki2*exp( - d*ki2))/(ki2 - km) -(( - km*exp( - d*km) + ki1*exp( - d*ki1)))/(ki1 - km )
}

FUNCTION search() (ms) {
	LOCAL x, t1, t2
	
	
    
    
    
    
    
    
    
	t1=1
	flag=0
    if(deriv(t1) < 0  ){
        while ( deriv(t1) < 0 && t1>1e-9 ){
            t2=t1
            t1=t1/10
        }
        if (deriv(t1) < 0) { 
    		printf("Error wrong deriv(t1)
    		flag=1
            search = 1e-9
    	}
    }else{
        t2=t1
        while (deriv(t2) > 0 && t2<1e9 ){
            t1=t2
            t2=t2*10
        }
    	if (deriv(t2) > 0) { 
    		printf("Error wrong deriv(t2)
    		flag=1
    		search = 1e9
    	}
    }
    
    
    
    
    while (t2 - t1 > 1e-6 && flag==0) {
		search = (t1+t2)/2
		x = deriv(search)
		if (x > 0) {
			t1 = search
		}else{
			t2 = search
		}
        
	}			
}


INITIAL {
    
	fixprecondition()
	factors()
	e = 0
	i1 = 0
	i2 = 0
	m = 0
	t0 = t
	net_send(firetimebound(), 1)

	nself=0
	nexcite=0
	ninhibit=0
}

NET_RECEIVE (w) {
	newstates(t-t0)
	update()	
    
	if (m > 1-eps) { 
        
		net_event(t)
		m = 0
	}
	if (flag == 1) { 
		nself = nself + 1
		net_send(firetimebound(), 1)
	}else{
		if (w > 0) {
			nexcite = nexcite + 1
			e = e + w
		}else{
			ninhibit = ninhibit + 1
			i1 = i1 + w
		}
        
		net_move(firetimebound() + t)
	}
}

FUNCTION firetimebound() (ms) {
	LOCAL slope
	slope = -km*m + ae*e + ai2*i2
	if (slope <= 1e-9) {
		firetimebound = 1e9
	}else{
		firetimebound = (1 - m)/slope
	}
}