INDEPENDENT {t FROM 0 TO 1 WITH 1 (ms)}

NEURON {
    SUFFIX nothing
}

VERBATIM
#include "misc.h"

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
    

    char *filename;

    filename = gargstr(1);

    if (*filename && !access(filename, F_OK)) {
        _lfile_exist = 1;

    } else {
        
        errno = 0;

        _lfile_exist = 0;
    }
ENDVERBATIM
}

FUNCTION istmpobj () {
VERBATIM
  _listmpobj=hoc_is_tempobj_arg(1);
ENDVERBATIM  
}


FUNCTION sassign() {
VERBATIM
    FILE *pipein;
    char string[BUFSIZ], **strname, *syscall;

    strname = hoc_pgargstr(1);
    syscall = gargstr(2);

    if( !(pipein = popen(syscall, "r"))) {
        fprintf(stderr,"System call failed\n");
        return 0; 
    }
    
    if (fgets(string,BUFSIZ,pipein) == NULL) {
        fprintf(stderr,"System call did not return a string\n");
        pclose(pipein); return 0;
    }

    
    hoc_assign_str(strname, string);

    pclose(pipein);
    errno = 0;
    return 0;
ENDVERBATIM
}


FUNCTION dassign() {
VERBATIM
    FILE *pipein, *outfile;
    char *strname, *syscall;
    double num;

    strname = gargstr(1);
    syscall = gargstr(2);

    if ( !(outfile = fopen("dassign","w"))) {
        fprintf(stderr,"Can't open output file dassign\n");
        return 0; 
    }

    if( !(pipein = popen(syscall, "r"))) {
        fprintf(stderr,"System call failed\n");
        fclose(outfile); return 0; 
    }
    
    if (fscanf(pipein,"%lf",&num) != 1) {
        fprintf(stderr,"System call did not return a number\n");
        fclose(outfile); pclose(pipein); return 0; 
    }

    fprintf(outfile,"%s=%g\n",strname,num);
    fprintf(outfile,"system(\"rm dassign\")\n");

    fclose(outfile); pclose(pipein);
    errno = 0;
    return 0;
ENDVERBATIM
}



PROCEDURE nokill() {
VERBATIM
  signal(SIGHUP, SIG_IGN);
ENDVERBATIM
}


FUNCTION prtime () {
VERBATIM
_lprtime = clock();
ENDVERBATIM
}


FUNCTION now () {
VERBATIM
  _lnow = time((time_t*)0);
  _lnow -= (12784) * 24*60*60; 
ENDVERBATIM
}



FUNCTION sleepfor (sec) {
VERBATIM
  struct timespec ts;
  ts.tv_sec = (time_t)_lsec;
  ts.tv_nsec = ifarg(2)?(long)*getarg(2):(long)0;
  return (double) nanosleep(&ts,(struct timespec*)0);
ENDVERBATIM
}


PROCEDURE spitchar(c) {
VERBATIM
{	
  printf("%c", (int)_lc);
}
ENDVERBATIM
}


VERBATIM
static char *pmlc;
ENDVERBATIM

PROCEDURE mymalloc(sz) {
VERBATIM
{ 
  size_t x,y;
  x=(size_t)_lsz;
  pmlc=(char *)malloc(x);
  printf("Did %ld: %x\n",x,pmlc);
  y=(unsigned int)_lsz-1;
  pmlc[y]=(char)97;
  printf("WRITE/READ 'a': "); 
  printf("%c\n",pmlc[y]);
  if (ifarg(2)) free(pmlc); else printf("Use unmalloc() to free memory\n");
}
ENDVERBATIM
}

PROCEDURE unmalloc() {
VERBATIM
  free(pmlc);
ENDVERBATIM
}


FUNCTION hocgetc() {
VERBATIM
{	
  FILE* f = hoc_obj_file_arg(1);
  _lhocgetc = (double)getc(f);
}
ENDVERBATIM
}

PROCEDURE pwd() {
  VERBATIM
  {char cwd[1000],cmd[1200];
  getcwd(cwd, 1000);
  sprintf(cmd, "execute1(\"strdef cwd\")\n");         hoc_oc(cmd);
  sprintf(cmd, "execute1(\"cwd=\\\"%s\\\"\")\n",cwd); hoc_oc(cmd);
  }
  ENDVERBATIM
}