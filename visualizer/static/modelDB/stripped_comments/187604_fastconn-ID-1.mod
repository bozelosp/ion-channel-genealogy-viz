VERBATIM
# include <stdlib.h>
# include <stdio.h>
# include <math.h>
ENDVERBATIM

NEURON {
	SUFFIX nothing
}

VERBATIM
#ifndef NRN_VERSION_GTEQ_8_2_0
#include <stdint.h>
extern double* vector_vec();
extern int vector_capacity();
extern void* vector_arg();
extern uint32_t nrnRan4int(uint32_t* idx1, uint32_t idx2);
#endif
extern double get_x_pos(double gid, double gmin, double BinNumX, double BinNumYZ, double binSizeX);
extern double get_y_pos(double gid, double gmin, double BinNumY, double BinNumZ, double binSizeY);
extern double get_z_pos(double gid, double gmin, double BinNumZ, double binSizeZ, double ZHeight);

ENDVERBATIM

VERBATIM

static  double fastconn (void* vv) {
  int finalconn, ny, nz, nhigh, num_pre, num_post, gmin, gmax, steps, myflaggy, myi, postgmin, stepover;
  double *x, *y, *z, *high, a, b, c, nconv, ncell, axonal_extent;

	
	finalconn = vector_instance_px(vv, &x); 
											
											

	ny = vector_arg_px(1, &y); 
	nz = vector_arg_px(2, &z); 
	nhigh = vector_arg_px(3, &high); 
	
	
	gmin = y[0];	
	gmax = y[1];	
	num_pre = gmax - gmin + 1;	
	
	nconv = y[2];	
	ncell = y[3];	
	num_post = y[4];	
	axonal_extent = y[5];	
	steps = y[6];	
	a = y[7];		
	b = y[8];		
	c = y[9];		
	postgmin = y[24];	
	stepover = y[26];	

	myi=2+num_post;	
			
			

	

	
	
	

	

	double *prepos;    
	prepos = (double *)malloc(sizeof(double)*num_pre*3);

	double *postpos;    
	postpos = (double *)malloc(sizeof(double)*num_post*3);

	
	
	

	
	
	
	
	

	int cell;

	for (cell=0; cell<num_pre; cell++) {
		prepos[cell*3 + 0] = get_x_pos(cell+gmin, gmin, y[10], y[11]*y[12], y[13]);
		prepos[cell*3 + 1] = get_y_pos(cell+gmin, gmin, y[11], y[12], y[14]);
		prepos[cell*3 + 2] = get_z_pos(cell+gmin, gmin, y[12], y[15], y[16]);
	}

	for (cell=0; cell<num_post; cell++) {
		postpos[cell*3 + 0] = get_x_pos(z[cell], postgmin, y[17], y[18]*y[19], y[20]);
		postpos[cell*3 + 1] = get_y_pos(z[cell], postgmin, y[18], y[19], y[21]);
		postpos[cell*3 + 2] = get_z_pos(z[cell], postgmin, y[19], y[22], y[23]);
	}

	   
	double current_distance [steps], connection_distribution [steps], distribution_denominator, conndist;
	int step, feasible_conns_this_step [steps], desired_conns_this_step [steps];

	distribution_denominator = 0.0;
	int max_fraction_step; max_fraction_step=0;
	for (step=0; step<steps; step++) {
		current_distance[step] = axonal_extent*1.0*(step+1)/(steps); 
		
		connection_distribution[step] = exp(-a*((current_distance[step]-b)*1.0/c)*((current_distance[step]-b)*1.0/c));
		if (connection_distribution[step]>connection_distribution[max_fraction_step]) {
			max_fraction_step=step;
		}
		distribution_denominator = distribution_denominator + connection_distribution[step];
	}

	if (connection_distribution[max_fraction_step]/distribution_denominator*nconv < 0.5) { 
		for (step=0; step<steps; step++) {
			desired_conns_this_step[step] = round((2.0*connection_distribution[step]/distribution_denominator)*(nconv));
															
															
		}
	} else {
		for (step=0; step<steps; step++) {
			desired_conns_this_step[step] = round((connection_distribution[step]/distribution_denominator)*(nconv));
															
															
		}
	}

	   
	int m, n, i, q, goupto, rem, extra, szr, szp [steps];
	double distance_between;
	u_int32_t idx1, idx2, maxidx1; 
	maxidx1 = y[25];

	for (n=0; n<num_post; n++) { 
		int myx = (int)z[n]; 
		idx1 = high[n];	
						
						
						
						
						
						
		idx2 = myx;		


		double *sortedpos;    
		sortedpos = (double *)malloc(sizeof(double)*num_pre*steps);
		
		
		for (step=0; step< steps; step++) {
			szp [step]=0; 	
							
							
							
							
			feasible_conns_this_step[step] = desired_conns_this_step[step];		
		}
		
		double dist;
		for(m=0; m<num_pre; m++) { 
			
			distance_between = sqrt((1.0*prepos[m*3 +0] - postpos[n*3 +0])*(prepos[m*3 +0] - postpos[n*3 +0])+(prepos[m*3 +1] - postpos[n*3 +1])*(prepos[m*3 +1] - postpos[n*3 +1])+(prepos[m*3 +2] - postpos[n*3 +2])*(prepos[m*3 +2] - postpos[n*3 +2]));
			for (step=0; step< steps; step++) {
				

				if (distance_between<= current_distance[step]) 
				{
					sortedpos[szp [step]*steps + step] = m;
					
					szp [step]++;
					break;
				}
			}
		}

		

		
		
		
		
		
			
		rem=0;extra=0;
		for (step=0; step<steps; step++) {	
			szr = szp [step]; 
			if (feasible_conns_this_step[step] + rem> szr) { 
				rem=feasible_conns_this_step[step]+rem-szr;
				
				if (step>0) {
					for(i=1; i<=step; i++) {
						if (szp [step-i] > feasible_conns_this_step[step-i]) {
							if (szp [step-i] - feasible_conns_this_step[step-i]>rem) {
								extra = rem;
							} else {
								extra = szp [step-i] - feasible_conns_this_step[step-i];
							}
							feasible_conns_this_step[step-i] = feasible_conns_this_step[step-i] + extra;
							feasible_conns_this_step[step] = feasible_conns_this_step[step] - extra;
							rem = rem - extra;				
						}
					}
				}
				if (rem>0 && step<steps-1) { 
					for(i=step+1; i<steps; i++) {				
						if (szp [i] > feasible_conns_this_step[i]) {
							if (szp [i] - feasible_conns_this_step[i]>rem) {
								extra = rem;
							} else {
								extra = szp [i] - feasible_conns_this_step[i];
							}
							feasible_conns_this_step[i] = feasible_conns_this_step[i] + extra;
							feasible_conns_this_step[step] = feasible_conns_this_step[step] - extra;
							rem = rem - extra;
						}
					}
				}
			}
		}

		
	
		rem=0;
		for (step=0; step<steps; step++) {	
			if (feasible_conns_this_step[step]>0) { 
				
				
				
				
				szr = szp [step]; 
				int r[szr]; 
				for (i=0; i< szr; i++) { 
					r[i] =  sortedpos[i*steps + step];
					
				}

				
				int tmp;
				u_int32_t randi;
				for (i=0; i<szr-1; i++) {
					randi =  nrnRan4int(&idx1, idx2) % (u_int32_t)szr; 
					tmp = r [i];	
					r[i] = r[randi];
					r[randi] = tmp;
				}

				if (feasible_conns_this_step[step]>szr) { 	
										
					goupto=szr;			
				} else {
					goupto=feasible_conns_this_step[step];	
				}

				for (q=0; q<goupto; q++) { 	
											
											
					x [myi] = (r[q]+gmin)*1.0;				
					if (num_post>1) {
						x [myi+1*stepover] = (z[n])*1.0;	
						
					}
					myi++;
				}
			} 
		}
		x [2+n] = idx1 + 1;
		
		free(sortedpos);
	}
	x [0] = myi-2-num_post;	
					
					
					
	x [1] = num_post; 
	
	free(prepos);
	free(postpos);
	return finalconn;
}
ENDVERBATIM




PROCEDURE install_fastconn () {
	VERBATIM
	install_vector_method("fastconn", fastconn);
	ENDVERBATIM
}