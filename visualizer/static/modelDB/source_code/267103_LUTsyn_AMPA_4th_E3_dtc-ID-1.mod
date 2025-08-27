
COMMENT
LUTsyn_AMPA_4th_E3_dtc.mod

This file implements the LUTsyn synapse model described in (Pham, 2021) for AMPA receptors

August 15, 2021
Duy-Tan Jonathan Pham

This software is Copyright Â© 2021 The University of Southern
California. All Rights Reserved.

Permission to use, copy, modify, and distribute this software
and its documentation for educational, research and non-profit
purposes, without fee, and without a written agreement is
hereby granted, provided that the above copyright notice, this
paragraph and the following three paragraphs appear in all copies.

Permission to make commercial use of this software may be obtained by contacting:
USC Stevens Center for Innovation
University of Southern California
1150 S. Olive Street, Suite 2300
Los Angeles, CA 90115, USA

This software program and documentation are copyrighted by The
University of Southern California. The software program and
documentation are supplied "as is", without any accompanying
services from USC. USC does not warrant that the operation of the
program will be uninterrupted or error-free. The end-user understands
that the program was developed for research purposes and is advised
not to rely exclusively on the program for any reason.

IN NO EVENT SHALL THE UNIVERSITY OF SOUTHERN CALIFORNIA BE
LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT
OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE
UNIVERSITY OF SOUTHERN CALIFORNIA HAS BEEN ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE. THE UNIVERSITY OF SOUTHERN CALIFORNIA
SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE. THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS"
BASIS, AND THE UNIVERSITY OF SOUTHERN CALIFORNIA HAS NO OBLIGATIONS
TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
ENDCOMMENT


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
	i (nA)			: EPSC
	v1 (mV)			: membrane voltage (so that python can "access" these values)
	gain_array		: pointer to a 3D numpy array that's been flattened to a 1D array -- the Look Up Table itself
	basis_gain 		: aka first-order amplitude
	gain			: dynamically changing amplitude value based on the LUT
	tau1			: inter-pulse intervals (IPIs)
	tau2
	tau3
	t3				: timings of the last 4 pulses
	t2
	t1
	t0
	MEM 		: memory window (normalized by granularity aka ratio R)
	order       : indicates the "order" of the present pulse based on timing events
	index       : holds the one-dimensional index to access the LUT array 
	g			: conductance
	factor		: scaling factor
	wf
	wtc3		: weighting factor for triple exponential
	o1_spks		: counts the number of spikes (pulses) for every order category
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
	
	MEM =  300/gran      :3000/10 (memory window size divided by granularity)   
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
	
	: triple exponential model
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
	SOLVE update_taus		: update the IPI values
	SOLVE find_gain			: find the LUT amplitude value based on the present IPI values
	wf = weight*factor*gain 
	C = C + wf
	B = B + wf
	E = E + wf
}


PROCEDURE update_taus(){
	: record the last 4 input pulse times
	t3 = t2
	t2 = t1
	t1 = t0
	t0 = t
	
	: adding 0.5 changes the floor function to be a round() function
	tau1 = floor((t0 - t1)/gran + 0.5)
	tau2 = floor((t0 - t2)/gran + 0.5)
	tau3 = floor((t0 - t3)/gran + 0.5)
}


PROCEDURE find_gain(){

	: find what ORDER the current pulse is
	if ((t0 != 0) && (t1 != 0) && (t2 != 0) && (t3 != 0)) : fourth+ spike case
	{
		if (tau1 > (MEM - 3)) : first order case
		{
			order = 1
		}
		
		else if (tau2 > (MEM - 2)) : second order case
		{
			order = 2
		}
		
		else if (tau3 > (MEM - 1)) : third order case
		{
			order = 3
		}
		
		else : fourth order case
		{

			order = 4
		}
	}
	
	else if (t1 == 0)         : first spike case
	{
		order = 1
	}
	
	else if (t2 == 0)        : second spike case
	{
		if (tau1 > (MEM - 3)) : first order case
		{
			order = 1
		}
		
		else
		{
			order = 2
		}
	}
	
	else if (t3 == 0)        : third spike case
	{
		if (tau1 > (MEM - 3)) : first order case
		{
			order = 1
		}
		
		else if (tau2 > (MEM - 2)) : second order case
		{

			order = 2
		}
		
		else
		{

			order = 3
		}
	}


	: based on the order, use the appropriate time constants and normalization factor
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
	
	: convert the multi-dimensional indices into a single one-dimensional index
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
		
	: printf("Look-Up Index = %1.0f, gran = %1.1f, tau1 = %1.0f, tau2 = %1.0f, tau3 = %1.0f, t0 = %1.0f, t1 = %1.0f, t2 = %1.0f, t3 = %1.0f", index, gran, tau1, tau2, tau3,t0,t1,t2,t3)
		
	: use the index to access the look-up table
	 VERBATIM
	 gain = ((double*)_p_gain_array)[(int) index];
	 ENDVERBATIM
	 
	}
	
}


: convert multiple indices into one flattened index
FUNCTION convert_index (ind1, ind2, ind3){
	UNITSOFF
    convert_index = (pow(ind3,3) - 3*pow(ind3,2) + 2*(ind3))/6 + (pow(ind2,2) - (ind2))/2 + ind1 
	UNITSON
}


