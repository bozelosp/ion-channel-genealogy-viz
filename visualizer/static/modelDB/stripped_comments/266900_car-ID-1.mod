NEURON {
	SUFFIX car
	POINTER stim_i
	USEION ca READ cai, cao WRITE ica

        RANGE flag, curr, gcabar, m, h,ica,sh,count,delta2,vrun2,stim_moltCa
	RANGE  curr, inf, fac, tau
	
}

UNITS {
	(mA) = (milliamp)
	(mV) = (millivolt)
	(molar) = (1/liter)
	(mM) =	(millimolar)
	FARADAY = (faraday) (coulomb)
	R = (k-mole) (joule/degC)
}


ASSIGNED {               
	ica (mA/cm2)

        inf[2]
	tau[2]		(ms)
        v               (mV)
        celsius 	(degC)
		 
	cai          (mM)      
	cao             (mM)      
	stim_i
}


PARAMETER { 

         curr              
        gcabar = 0      (mho/cm2) 
        eca = 140      (mV)      
        eca2=140		(mV)      
        sh=0 (mV)     
	time1=600
	time0=100
	alpha=1
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

       
}  

STATE {	
	m 
	h 
}            


INITIAL {
	rates(v,sh2)
        m = 0    
	h = 1    
	vrun=0
	vvrun=vrun
}

BREAKPOINT {
	SOLVE states METHOD cnexp
	
	 eca2=eca+vvrun*alpha
	ica = gcabar*m*m*m*h*(v - eca)

}


DERIVATIVE states {
	rates(v,sh2)
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
   
PROCEDURE rates(v(mV),sh2) {LOCAL a, b ,i

      
		   
	FROM i=0 TO 1 {
		tau[i] = vartau(v,i)
		inf[i] = varss(v-sh2,i)
	}
}




FUNCTION varss(v(mV), i) {
	if (i==0) {
	    varss = 1 / (1 + exp((v+48.5)/(-3(mV)))) 
	}
	else if (i==1) {
             varss = 1/ (1 + exp((v+53)/(1(mV))))    
	}
}

FUNCTION vartau(v(mV), i) (ms){
	if (i==0) {
         vartau = 50(ms)  
   
        }
	else if (i==1) {
          vartau = 5(ms)   
     
       }
	
}