NEURON {
	POINT_PROCESS LUTsyn_AMPA_4th_E3_dtc
	RANGE gain, basis_gain, tau1, tau2, tau3, t3, t2, t1, t0, MEM, scalar, order, open, index, gran
	RANGE g, e, C, B, E, factor, wf, v1
	RANGE tc1, tc2, tc3, wtc2, wtc3
	RANGE o1_tc1, o1_tc2, o1_tc3, o2_tc1, o2_tc2, o2_tc3, o3_tc1, o3_tc2, o3_tc3, o4_tc1, o4_tc2, o4_tc3
	RANGE o1_spks, o2_spks, o3_spks, o4_spks
	RANGE factor1, factor2, factor3, factor4
	POINTER gain_array
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(nS) = (nanosiemens)
}

ASSIGNED {
	i (nA)			
	v1 (mV)			
	gain_array		
	basis_gain 		
	gain			
	tau1			
	tau2
	tau3
	t3				
	t2
	t1
	t0
	MEM 		
	order       
	index       
	g			
	factor		
	wf
	wtc3		
	o1_spks		
	o2_spks
	o3_spks
	o4_spks
}

PARAMETER {
	scalar = 1
	gran = 1
	v
	e = 0
	tc1 = 0.86 (ms)
	tc2 = 4.24 (ms)
	tc3 = 12 (ms)
	wtc2 = 0.5
	
	o1_tc1 = 1 (ms)
	o1_tc2 = 1 (ms)
	o1_tc3 = 1 (ms)
	factor1 = 1
	
	o2_tc1 = 1 (ms)
	o2_tc2 = 1 (ms)
	o2_tc3 = 1 (ms)
	factor2 = 1

	o3_tc1 = 1 (ms)
	o3_tc2 = 1 (ms)
	o3_tc3 = 1 (ms)
	factor3 = 1
	
	o4_tc1 = 1 (ms)
	o4_tc2 = 1 (ms)
	o4_tc3 = 1 (ms)
	factor4 = 1

}

INITIAL {
	if (tc1/tc2 > .9999) {
		tc1 = .9999*tc2
	}
	if (tc2/tc3 > .9999) {
	    tc2 = .9999*tc3
    }
	
    wtc3 = 1-wtc2
	
	C = 0
	B = 0
	E = 0
	
	tau1 = -1
	tau2 = -1
	tau3 = -1

	t3 = 0
	t2 = 0
	t1 = 0
	t0 = 0
	
	o1_spks = 0
	o2_spks = 0
	o3_spks = 0
	o4_spks = 0
	
	MEM =  300/gran      
	order = 0
	index = 0        
	basis_gain = basis_gain / fabs(scalar)
	gain = basis_gain
}

STATE{
	C
	B
	E
}

BREAKPOINT {
	SOLVE state METHOD cnexp	
	
	
	g = wtc2*B + wtc3*E - C  
	i = g * (v-e)
	v1 = v
}

DERIVATIVE state {
	C' = -C/tc1
	B' = -B/tc2
	E' = -E/tc3
}

NET_RECEIVE(weight (uS)) {
	SOLVE update_taus		
	SOLVE find_gain			
	wf = weight*factor*gain 
	C = C + wf
	B = B + wf
	E = E + wf
}


PROCEDURE update_taus(){
	
	t3 = t2
	t2 = t1
	t1 = t0
	t0 = t
	
	
	tau1 = floor((t0 - t1)/gran + 0.5)
	tau2 = floor((t0 - t2)/gran + 0.5)
	tau3 = floor((t0 - t3)/gran + 0.5)
}


PROCEDURE find_gain(){

	
	if ((t0 != 0) && (t1 != 0) && (t2 != 0) && (t3 != 0)) 
	{
		if (tau1 > (MEM - 3)) 
		{
			order = 1
		}
		
		else if (tau2 > (MEM - 2)) 
		{
			order = 2
		}
		
		else if (tau3 > (MEM - 1)) 
		{
			order = 3
		}
		
		else 
		{

			order = 4
		}
	}
	
	else if (t1 == 0)         
	{
		order = 1
	}
	
	else if (t2 == 0)        
	{
		if (tau1 > (MEM - 3)) 
		{
			order = 1
		}
		
		else
		{
			order = 2
		}
	}
	
	else if (t3 == 0)        
	{
		if (tau1 > (MEM - 3)) 
		{
			order = 1
		}
		
		else if (tau2 > (MEM - 2)) 
		{

			order = 2
		}
		
		else
		{

			order = 3
		}
	}


	
	if (order == 1)
	{
	 gain = basis_gain
	 tc1 = o1_tc1
	 tc2 = o1_tc2
	 tc3 = o1_tc3
	 factor = factor1
	 o1_spks = o1_spks + 1
	}
	
	else 
	
	
	{
		if (order == 2)
		{
		 index = convert_index(tau1,MEM-2,MEM-1)
		 tc1 = o2_tc1
		 tc2 = o2_tc2
		 tc3 = o2_tc3
		 factor = factor2
		 o2_spks = o2_spks + 1
		}
		
		else if (order == 3)
		{
		 index = convert_index(tau1,tau2,MEM-1)
		 tc1 = o3_tc1
		 tc2 = o3_tc2
		 tc3 = o3_tc3
		 factor = factor3
		 o3_spks = o3_spks + 1
		}
		
		else if (order == 4)
		{
		 index = convert_index(tau1,tau2,tau3)
		 tc1 = o4_tc1
		 tc2 = o4_tc2
		 tc3 = o4_tc3
		 factor = factor4
		 o4_spks = o4_spks + 1
		}
		
	
		
	
	 VERBATIM
	 gain = ((double*)_p_gain_array)[(int) index];
	 ENDVERBATIM
	 
	}
	
}



FUNCTION convert_index (ind1, ind2, ind3){
	UNITSOFF
    convert_index = (pow(ind3,3) - 3*pow(ind3,2) + 2*(ind3))/6 + (pow(ind2,2) - (ind2))/2 + ind1 
	UNITSON
}