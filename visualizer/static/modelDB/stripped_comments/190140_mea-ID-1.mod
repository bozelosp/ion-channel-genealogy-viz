NEURON {
	SUFFIX mea
	POINTER transmembrane_current_m
	RANGE mea_line0,mea_line1,mea_line2,mea_line3,mea_line4,mea_line5,mea_line6,mea_line7,mea_line8,mea_line9,mea_line10,mea_line11,mea_line12,mea_line13,mea_line14,mea_line15
	RANGE initial_part_line0,initial_part_line1,initial_part_line2,initial_part_line3,initial_part_line4,initial_part_line5,initial_part_line6,initial_part_line7,initial_part_line8,initial_part_line9,initial_part_line10,initial_part_line11,initial_part_line12,initial_part_line13,initial_part_line14,initial_part_line15	

}

PARAMETER {
	

	}

ASSIGNED {
 
	transmembrane_current_m 
	initial_part_line0
	initial_part_line1
	initial_part_line2
	initial_part_line3
	initial_part_line4
	initial_part_line5
	initial_part_line6
	initial_part_line7
	initial_part_line8
	initial_part_line9
	initial_part_line10
	initial_part_line11
	initial_part_line12
	initial_part_line13
	initial_part_line14
	initial_part_line15

	mea_line0
	mea_line1
	mea_line2
	mea_line3
	mea_line4
	mea_line5
	mea_line6
	mea_line7
	mea_line8
	mea_line9
	mea_line10
	mea_line11
	mea_line12
	mea_line13
	mea_line14
	mea_line15


}

BREAKPOINT { 

	
	mea_line0 = transmembrane_current_m * initial_part_line0 * 1e-1 	
	mea_line1 = transmembrane_current_m * initial_part_line1 * 1e-1 
	mea_line2 = transmembrane_current_m * initial_part_line2 * 1e-1 
	mea_line3 = transmembrane_current_m * initial_part_line3 * 1e-1 
	mea_line4 = transmembrane_current_m * initial_part_line4 * 1e-1 
	mea_line5 = transmembrane_current_m * initial_part_line5 * 1e-1 
	mea_line6 = transmembrane_current_m * initial_part_line6 * 1e-1 
	mea_line7 = transmembrane_current_m * initial_part_line7 * 1e-1 
	mea_line8 = transmembrane_current_m * initial_part_line8 * 1e-1 
	mea_line9 = transmembrane_current_m * initial_part_line9 * 1e-1 
	mea_line10 = transmembrane_current_m * initial_part_line10 * 1e-1 
	mea_line11 = transmembrane_current_m * initial_part_line11 * 1e-1 
	mea_line12 = transmembrane_current_m * initial_part_line12 * 1e-1 
	mea_line13 = transmembrane_current_m * initial_part_line13 * 1e-1 
	mea_line14 = transmembrane_current_m * initial_part_line14 * 1e-1 
	mea_line15 = transmembrane_current_m * initial_part_line15 * 1e-1 

}