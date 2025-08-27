UNITS {
        (S) = (siemens)
        (mA) = (milliamp)
        (mV) = (millivolt)
}
 
NEURON {
        SUFFIX leak
        NONSPECIFIC_CURRENT ileak 
        RANGE gleakbar,ileak,eleak

}
 
PARAMETER {
        v   (mV)
        dt  (ms)
		gleakbar  = 0.00002  (S/cm2)
        eleak = -55 (mV)
		
}
 
 
ASSIGNED {
        ileak (mA/cm2)
        
}
 
BREAKPOINT {
         ileak = gleakbar*(v - eleak)  

		
}