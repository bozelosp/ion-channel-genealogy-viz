INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX nothing
}

VERBATIM
#include <unistd.h>     
#include <errno.h>      
#include <signal.h>
#include <sys/types.h>         
#include <time.h>
#include <stdio.h>
#include <limits.h>
ENDVERBATIM


FUNCTION file_exist() {
VERBATIM
    

    char *gargstr(), *filename;

    filename = gargstr(1);

    if (*filename && !access(filename, F_OK)) {
        _lfile_exist = 1;

    } else {
        
        errno = 0;

        _lfile_exist = 0;
    }
ENDVERBATIM
}


PROCEDURE sassign() {
VERBATIM
    FILE *pipein;
    char string[BUFSIZ], **strname, *syscall;
    char** hoc_pgargstr();

    strname = hoc_pgargstr(1);
    syscall = gargstr(2);

    if( !(pipein = popen(syscall, "r"))) {
        fprintf(stderr,"System call failed\n");
        return; 
    }
    
    if (fgets(string,BUFSIZ,pipein) == NULL) {
        fprintf(stderr,"System call did not return a string\n");
        pclose(pipein); return;
    }

    
    hoc_assign_str(strname, string);

    pclose(pipein);
    errno = 0;
ENDVERBATIM
}


PROCEDURE dassign() {
VERBATIM
    FILE *pipein, *outfile;
    char *strname, *syscall;
    double num;

    strname = gargstr(1);
    syscall = gargstr(2);

    if ( !(outfile = fopen("dassign","w"))) {
        fprintf(stderr,"Can't open output file dassign\n");
        return; 
    }

    if( !(pipein = popen(syscall, "r"))) {
        fprintf(stderr,"System call failed\n");
        fclose(outfile); return; 
    }
    
    if (fscanf(pipein,"%lf",&num) != 1) {
        fprintf(stderr,"System call did not return a number\n");
        fclose(outfile); pclose(pipein); return; 
    }

    fprintf(outfile,"%s=%g\n",strname,num);
    fprintf(outfile,"system(\"rm dassign\")\n");

    fclose(outfile); pclose(pipein);
    errno = 0;
ENDVERBATIM
}



PROCEDURE nokill() {
VERBATIM
  signal(SIGHUP, SIG_IGN);
ENDVERBATIM
}


FUNCTION prtime () {
VERBATIM
  double prt;
  static double PRTIME;
  prt = (clock()-PRTIME)/CLOCKS_PER_SEC;
  
  if (prt<0) prt += UINT_MAX/CLOCKS_PER_SEC; 
  PRTIME=clock();
  _lprtime = prt;
ENDVERBATIM
}


FUNCTION now () {
VERBATIM
  _lnow = time((time_t*)0);
  _lnow -= (12784) * 24*60*60; 
ENDVERBATIM
}


PROCEDURE sleepfor (sec) {
VERBATIM
  struct timespec ts;
  ts.tv_sec = (time_t)_lsec;
  ts.tv_nsec = (long)0;
  nanosleep(&ts,(struct timespec*)0);
ENDVERBATIM
}


PROCEDURE fspitchar(c) {
VERBATIM
{	
  FILE* f, *hoc_obj_file_arg();
  f = hoc_obj_file_arg(2);
  fprintf(f, "%c", (int)_lc);
}
ENDVERBATIM
}


PROCEDURE spitchar(c) {
VERBATIM
{	
  printf("%c", (int)_lc);
}
ENDVERBATIM
}


FUNCTION hocgetc() {
VERBATIM
{	
  FILE* f, *hoc_obj_file_arg();
  f = hoc_obj_file_arg(1);
  _lhocgetc = (double)getc(f);
}
ENDVERBATIM
}