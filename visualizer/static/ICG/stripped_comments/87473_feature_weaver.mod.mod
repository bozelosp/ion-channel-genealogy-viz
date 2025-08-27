NEURON {
	SUFFIX nothing
}

VERBATIM


extern double hoc_epsilon;
#define EPS hoc_epsilon
static double lin_interp();


static double ywnscl_fitness(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
	int val_pkind, pkind, mod_start, mod_end, mnind, mxind, n_pts;
	double ytmp, val_winsz, mod_winsz;

	
	
	ny = vector_instance_px(vv, &y);		
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	val_pkind = *getarg(5);
	pkind = *getarg(6);
	mnind = *getarg(7);
	mxind = *getarg(8);

	j = 0;
	sum = 0.;
	n_pts = 0;

	
	
	
	

        
        
        mod_start = mxind+1;
        mod_end = mnind-1;
        i = mnind;
        j = 0; 

	while( i < mxind && j < nyval ) {

	    
	    

	    
	    if( xval[j] > (x[i] - xpeak) ) {
	        while( (i<ny) && xval[j] > x[i] - xpeak && (fabs(xval[j] - (x[i] - xpeak)) > EPS) ) {
	            i++;
		}
	    }


	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        
	        d = y[i] - yval[j];
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
	    } else if((i == ny-1) && (j==nyval-1) && (xval[j] > x[i]-xpeak)) {
		
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
	    } else if((i==ny-1) && (j < nyval) ) {
		
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
	    } else {
	        
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
	        d = ytmp - yval[j];
		sum += d*d;
                if( n_pts == 0 ) { mod_start = i; }
		n_pts++;
	    }

	    j++;
	}

	
	mod_end = i;
	if( mod_end == nx ) { 
	     mod_end--; 
	}

	
        if( n_pts == 0 ) { 
	    sum = -1;
        } else {
  	    sum = sqrt(sum/(double)n_pts);
	}

        
	
	val_winsz = xval[nxval-1] - xval[0];
        if( n_pts == 0 ) { 
	    mod_winsz = 0;
	} else {
	    mod_winsz = x[mod_end] - x[mod_start];
	}
        if( mod_winsz > 0 ) {
	    sum = mod_winsz*sum/val_winsz; 
	} else {
	    sum = -1.0;
	}
	return sum;

}



static double yfitness_weaver(void* vv) {
	int nx, ny, nyval, nxval, i, j;
	double sum, d, xpeak, *y, *x, *yval, *xval;
	int val_pkind, pkind;
	double ytmp;
	ny = vector_instance_px(vv, &y);
	nx = vector_arg_px(1, &x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	xpeak = *getarg(2);
	nyval = vector_arg_px(3, &yval);
	nxval = vector_arg_px(4, &xval);
	val_pkind = *getarg(5);
	pkind = *getarg(6);
	j = 0;
	sum = 0.;

	
	
	
	
	
	

	
	d = y[pkind] - yval[val_pkind];
	sum += d*d;

	
        i = pkind; 
	j = val_pkind -1;
	while( i > 0 && j >= 0 ) {

	    
	    

	    
	    if( xval[j] < x[i] - xpeak ) {
	        while( (i>0) && xval[j] < x[i] - xpeak ) {

	            i--;
		}
	    }


	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        
	        d = y[i] - yval[j];
		sum += d*d;
	    } else if( j==0 && i==0 && (xval[j] < x[i]-xpeak) ) {
		
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
	    } else if( i==0 && j>0 ) {
		
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
	    } else {
	        
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i+1]-xpeak,y[i+1]);
	        d = ytmp - yval[j];
		sum += d*d;
	    }

	    j--;
	}

	
        i = pkind; 
	j = val_pkind + 1;
	while( i < ny-1 && j < nyval ) {

	    
	    

	    
	    if( xval[j] > x[i] - xpeak ) {
	        while( (i<ny) && xval[j] > x[i] - xpeak ) {
	            i++;
		}
	    }



	    if( fabs(xval[j] - (x[i]-xpeak)) < EPS ) {
	        
	        d = y[i] - yval[j];
		sum += d*d;
	    } else if((i == ny-1) && (j==nyval-1) && (xval[j] > x[i]-xpeak)) {
		
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j-1],yval[j-1]);
		d = y[i] - ytmp;
		sum += d*d;
	    } else if((i==ny-1) && (j < nyval) ) {
		
	        ytmp = lin_interp(x[i]-xpeak,xval[j],yval[j],xval[j+1],yval[j+1]);
		d = y[i] - ytmp;
		sum += d*d;
	    } else {
	        
	        ytmp = lin_interp(xval[j],x[i]-xpeak,y[i],x[i-1]-xpeak,y[i-1]);
	        d = ytmp - yval[j];
		sum += d*d;
	    }

	    j++;
	}

	return sum;

	
}




static double lin_interp(xstar, x1, y1, x2, y2) double xstar, x1, y1, x2, y2; {
        double ystar;

	if( fabs(x2-x1) < EPS ) return 0.5*(y1+y2);

	ystar = y1 + ((y2-y1)/(x2-x1))*(xstar-x1);
	return ystar;
}


static double firstmax(void* vv) { 
	int ny, i;
	double *y;
	ny = vector_instance_px(vv, &y) - 1;
	i = 0;
	while (i < ny) {
		if (y[i] > y[i+1]) {
			return (double) i;
		}
		i = i + 1;
	}
	return 0.;
}


static double nextpeak(void* vv) {
	int ny, i;
	double *y;
	ny = vector_instance_px(vv, &y) - 1;
	i = *getarg(1);
	while (i < ny) {
		if (y[i] >= -20) {
			if (y[i] > y[i+1]) {
				return (double) i;
			}
			i = i + 1;
		} else {
			i = i + 2;
		}
	}
	return 0.;
}
 

static	float	sqrarg;

#define	SQR(a) (sqrarg=(a),sqrarg*sqrarg)



void	linfit(void* vv) {

  
	float	*x, *y;
	int	strt, end;
	float	*a, *b, *siga, *sigb, *chi2;

	int	i,nx,ny;
	float	wt, t, sxoss, ss, sigdat;
	float	sx  = 0.0;
	float	sy  = 0.0;
	float	st2 = 0.0;

	ny = vector_instance_px(vv,&y);
	nx = vector_instance_px(1,&x);
	if (nx != ny) { hoc_execerror("vectors not same size", 0); }
	strt = *getarg(2);
	end  = *getarg(3);
	*a    = *getarg(4);
	*b    = *getarg(5);
	*siga = *getarg(6);
	*sigb = *getarg(7);
	*chi2 = *getarg(8);

	*b = 0.0;
	for( i = strt;  i <= end;  i++ )
	{
		sx += x[i];	sy += y[i];
	}
	ss = end+1-strt;

	sxoss = sx/ss;

	for( i = strt;  i <= end;  i++ )
	{
		t = x[i] - sxoss;
		st2 += t*t;
		*b += t*y[i];
	}

	*b /= st2;
	*a = (sy-sx*(*b))/ss;
	*siga = sqrt((1.0 + sx*sx/(ss*st2))/ss);
	*sigb = sqrt(1.0/st2);
	*chi2 = 0.0;

	for( i = strt;  i <= end;  i++ ) {
		*chi2 += SQR(y[i]-(*a)-(*b)*x[i]);
	}
	sigdat = sqrt((*chi2)/(end-2));
	*siga *= sigdat;
	*sigb *= sigdat;
}
ENDVERBATIM


PROCEDURE install_weaver_fitness() {
VERBATIM
  {static int once; if (!once) { once = 1;
	install_vector_method("ywnscl_fitness", ywnscl_fitness);
	install_vector_method("yfitness_weaver", yfitness_weaver);
	install_vector_method("firstmax", firstmax);
	install_vector_method("nextpeak", nextpeak);
	install_vector_method("linfit", linfit);
  }}
ENDVERBATIM
}