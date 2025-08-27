
COMMENT
LUTsyn_NMDA_5th_E3_dtc.mod

This file implements the LUTsyn synapse model described in (Pham, 2021) for NMDA receptors

July 13, 2021
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
	POINT_PROCESS LUTsyn_NMDA_5th_E3_dtc
	RANGE gain, basis_gain, tau1, tau2, tau3, tau4, t4, t3, t2, t1, t0, MEM, scalar, order, open, gran, index
	RANGE g0, g1, g2, alpha, g, e, open_Mg, nbNMDAR, C, B, E, factor, wf, v1
	RANGE tc1, tc2, tc3, wtc2, wtc3
	RANGE o1_tc1, o1_tc2, o1_tc3, o2_tc1, o2_tc2, o2_tc3, o3_tc1, o3_tc2, o3_tc3, o4_tc1, o4_tc2, o4_tc3, o5_tc1, o5_tc2, o5_tc3
	RANGE factor1, factor2, factor3, factor4, factor5
	RANGE o1_spks, o2_spks, o3_spks, o4_spks, o5_spks
	POINTER gain_array
	NONSPECIFIC_CURRENT i
}

UNITS {
	(nA) = (nanoamp)
	(mV) = (millivolt)
	(uS) = (microsiemens)
}

ASSIGNED {
	i (nA)         : EPSC
	gain_array     : pointer to look-up table
	basis_gain	   : aka first-order amplitude
	gain 		   : dynamically changing amplitude value based on the LUT
	tau1           : inter-pulse intervals (IPIs)
	tau2
	tau3
	tau4
	t4             : timings of the last 5 pulses
	t3
	t2
	t1
	t0
	order       : indicates the "order" of the present pulse based on timing events
	open		: open-state probability
	gran		: granularity in ms 
	MEM         : memory window (normalized by granularity aka ratio R)
	index       : holds the one-dimensional index to access the LUT array 
	g			: conductance
	factor		: scaling factor
	wf
	v1			: membrane voltage (so that python can "access" these values)
	wtc3		: weighting factor for triple exponential
	o1_spks		: counts the number of spikes (pulses) for every order category
	o2_spks
	o3_spks
	o4_spks
	o5_spks
}

PARAMETER {
	scalar = 1
	alpha = 0.01
	g0 = 0.0
	g1 = 40.0   (pS)
	g2 = 247.0  (pS)
	Mg = 1
	open_Mg = 0.0
	v
	e = 0
	nbNMDAR = 1.235
	tc1 = 0.1 (ms)
	tc2 = 10 (ms)
	tc3 = 11 (ms)
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
	
	o5_tc1 = 1 (ms)
	o5_tc2 = 1 (ms)
	o5_tc3 = 1 (ms)
	factor5 = 1
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
	tau4 = -1
	
	t4 = 0
	t3 = 0
	t2 = 0
	t1 = 0
	t0 = 0
	
	MEM =  200       :1000/5 (memory window size divided by granularity)   
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
	open = wtc2*B + wtc3*E - C  

	: simulate the magnesium blockade
	g0 = g1 + ((g2 - g1) / (1 + exp(alpha * v * 0.001(1/mV))))
	open_Mg = (open) / (1 + exp(-62 * v * 0.001(1/mV)) * Mg / 3.57)
	g = g0 * open_Mg * 0.001
	i = g * (v-e) * nbNMDAR
	v1 = v
}

DERIVATIVE state {
	C' = -C/tc1
	B' = -B/tc2
	E' = -E/tc3
}

NET_RECEIVE(weight (uS)) {
	SOLVE update_taus    : update the IPI values
	SOLVE find_gain      : find the LUT amplitude value based on the present IPI values
	wf = weight*factor*gain 
	C = C + wf
	B = B + wf
	E = E + wf
}


PROCEDURE update_taus(){
	: record the last 5 input pulse times
	t4 = t3
	t3 = t2
	t2 = t1
	t1 = t0
	t0 = t
	
	: adding 0.5 changes the floor function to be a round() function
	   : convert taus from time delays (ms) into array indices (modified by granularity)
	tau1 = floor((t0 - t1)/gran + 0.5)
	tau2 = floor((t0 - t2)/gran + 0.5)
	tau3 = floor((t0 - t3)/gran + 0.5)
	tau4 = floor((t0 - t4)/gran + 0.5)
}


PROCEDURE find_gain(){

	: find what ORDER the current pulse is
	if ((t0 != 0) && (t1 != 0) && (t2 != 0) && (t3 != 0) && (t4 != 0)) : fifth+ spike case
	{
		if (tau1 > (MEM - 4)) : first order case
		{
			order = 1
		}
		
		else if (tau2 > (MEM - 3)) : second order case
		{
			order = 2
		}
		
		else if (tau3 > (MEM - 2)) : third order case
		{
			order = 3
		}
		
		else if (tau4 > (MEM-1)): fourth order case
		{
			order = 4
		}
		
		else : fifth order+ case
		{
			order = 5
		}
	}
	
	else if (t1 == 0)         : first spike case
	{
		order = 1
	}
	
	else if (t2 == 0)        : second spike case
	{
		if (tau1 > (MEM - 4)) : first order case
		{
			order = 1
		}
		
		else : second order case 
		{
			order = 2
		}
	}
	
	else if (t3 == 0)        : third spike case
	{
		if (tau1 > (MEM - 4)) : first order case
		{
			order = 1
		}
		
		else if (tau2 > (MEM - 3)) : second order case
		{
			order = 2
		}
		
		else : third order case
		{
			order = 3
		}
	}

	else if (t4 == 0)        : fourth spike case
	{
		if (tau1 > (MEM - 4)) : first order case
		{
			order = 1
		}
		
		else if (tau2 > (MEM - 3)) : second order case
		{
			order = 2
		}
		
		else if (tau3 > (MEM - 2)) : third order case
		{
			order = 3
		}
		
		else :FOURTH ORDER CASE 
		{
			order = 4
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
		 index = convert_index(tau1,MEM-3,MEM-2,MEM-1)
		 tc1 = o2_tc1
		 tc2 = o2_tc2
		 tc3 = o2_tc3
		 factor = factor2
		 o2_spks = o2_spks + 1
		}
		
		else if (order == 3)
		{
		 index = convert_index(tau1,tau2,MEM-2,MEM-1)
		 tc1 = o3_tc1
		 tc2 = o3_tc2
		 tc3 = o3_tc3
		 factor = factor3
		 o3_spks = o3_spks + 1
		}
		
		else if (order == 4)
		{
		 index = convert_index(tau1,tau2,tau3,MEM-1)
		 tc1 = o4_tc1
		 tc2 = o4_tc2
		 tc3 = o4_tc3
		 factor = factor4
		 o4_spks = o4_spks + 1
		}
		
		else if (order == 5)
		{
		 index = convert_index(tau1,tau2,tau3,tau4)
		 tc1 = o5_tc1
		 tc2 = o5_tc2
		 tc3 = o5_tc3
		 factor = factor5
		 o5_spks = o5_spks + 1
		}
		
	: use the index to access the look-up table
	 VERBATIM
	 gain = ((double*)_p_gain_array)[(int) index];
	 ENDVERBATIM
 
	 
	}
}


: convert multiple indices into one flattened index
FUNCTION convert_index (ind1, ind2, ind3, ind4){
	UNITSOFF
    convert_index = ((pow(ind4,4) - 6*pow(ind4,3) + 11*pow(ind4,2) - 6*(ind4))/24 + (pow(ind3,3) - 3*pow(ind3,2) + 2*(ind3))/6 + (pow(ind2,2) - (ind2))/2 + ind1 )
	UNITSON
}


