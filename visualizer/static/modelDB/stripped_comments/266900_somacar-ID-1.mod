NEURON {
	SUFFIX somacar
	POINTER stim_i
	USEION ca READ eca WRITE ica
        RANGE flag, curr, gcabar, m, h,sh,count,delta2,vrun2
	RANGE flag,  inf, fac, tau,flag,eca2
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
}



PARAMETER { curr 
        sh = 0         (mV)      
        eca = 140       (mV)      
        gcabar = 0      (mho/cm2) 
        celsius =34        (degC)
		 count=1
	vrun (mV)
	delta=0
	vinit=-76.2
	alpha=1.1
	sh2
	alphash1=0.15
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
	eca2=140(mV)      
	
}


ASSIGNED {      
        v               (mV)
 	ica             (mA/cm2)
	ecar             (mV)      
        inf[2]
	fac[2]
	tau[2]
       
        
		
	stim_i
		
}

STATE {	
	m 
	h 
} 


INITIAL {
	m = 0    
	h = 1    
	vrun=0
	vvrun=vrun
	rates(v)
      ica = gcabar*m*m*m*h*(v - eca)  
}

BREAKPOINT {
	SOLVE states METHOD cnexp
        
    eca2=eca+vvrun*alpha	
	ica = gcabar*m*m*m*h*(v - eca)
}


DERIVATIVE states {
        rates(v)
	m' = (inf[0]-m)/tau[0]
	h' = (inf[1]-h)/tau[1]
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

PROCEDURE rates(v(mV)) { LOCAL ica,i

		   
	FROM i=0 TO 1 {
		tau[i] = vartau(i)
		inf[i] = varss(v-sh2,i)
	}
}



FUNCTION varss(v(mV), i) {
	if (i==0) {
	   varss = 1 / (1 + exp((v+60)/(-3))) 
	}
	else if (i==1) {
           varss = 1/ (1 + exp((v+62)/(1)))   
	}
}

FUNCTION vartau(i) {
	if (i==0) {
           vartau = 100  
        }
	else if (i==1) {
           vartau = 5    
       }
	
}