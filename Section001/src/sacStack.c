/*********************************************************************
*	sacStack.c:
*	 stack sac traces aligned with a time mark.
*
*	Author:  Lupei Zhu
*
*	Revision History
*		March 1999	Initial coding.
*		June  2000	include other time marks option.
*		Oct 26, 2000	add normalization option.
*		July 25, 2008	add unset event location in the sac head.
*********************************************************************/

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "sac_lpz.h"
#include "Complex.h"

int main(int argc, char **argv) {
  int 		i,j,nn,shift,error,start,end,ntrace,nrdc,norm,setbaz,n_files,isqr,wndw,setT,unseteve,intg, smth;
  char		outf[64], line[512], rdc, sac_file[128];
  float 	ar,dt,*src,*trace,reduce_vel,align,maxA,baz,p,t0,t1,t2;
  SACHEAD	hd,*hd0;
  
  error = 0;
  isqr = 0;
  wndw = 0;
  rdc = 't'; nrdc = -5;		/* default: align with begining */
  norm = 0;			/* no nomalization */
  outf[0] = '\0';
  setbaz = 0;
  setT = 0;
  unseteve = 0;
  intg = 0;
  smth = 0;
  n_files=0;
  /* input parameters */
  for (i=1; !error && i < argc; i++) {
    if (argv[i][0] == '-') {
       switch(argv[i][1]) {
  
       case 'E':
	 rdc=argv[i][2];
	 if (rdc == 't') sscanf(&argv[i][3],"%d",&nrdc);
	 else sscanf(&argv[i][2],"%f",&reduce_vel);
	 break;

       case 'I':
         intg = 1;
	 break;

       case 'M':
         sscanf(&argv[i][2],"%d",&smth);
	 break;

       case 'N':
	 norm = 1;
	 if (argv[i][2] == 'a') norm = 2;
	 break;
  
       case 'O':
	 strcpy(outf,&argv[i][2]);
         break;

       case 'Q':
         isqr = 1.;
         break;

       case 'R':
	 wndw = 1;
	 sscanf(&argv[i][2],"%f/%f",&t1,&t2);
	 break;

       case 'S':
	 setbaz = 1;
	 sscanf(&argv[i][2],"%f/%f",&baz,&p);
	 break;

       case 'T':
	 setT = 1;
	 break;

       case 'U':
	 unseteve = 1;
	 break;

       default:
	 error = 1;
	 break;

       }
    }
    else n_files++;
  }

  if (argc < 2 || error || outf[0] == '\0') {
     fprintf(stderr,"usage: %s [-E(t(0-9,-5(b),-3(o),-2(a))|vel)] [-I] [-N[a]] [-Mn] [-Q] [-Rt1/t2] [-Sbaz/p] [-T] [-U] -Ooutput_file (sac_traces in the argument list or from the stdin)\n\
		-E: align with a time mark or with an apparent velocity (b)\n\
		-I: integrate the trace before stacking (off)\n\
		-M: do n-point moving average of the trace before stacking (0)\n\
		-N: normalize by amplitude, append a to normalize by area (off)\n\
		-Q: square traces before stacking (off)\n\
		-R: time window t1 and t2\n\
		-S: set baz and user0 (p) (average)\n\
		-T: set the ref. time to the alignment time\n\
		-U: unset event location in the header (station location)\n",argv[0]);
     return -1;
  }

  ntrace = 0;
  i = 0;
  while( (n_files && ++i<argc) || (n_files==0 && fgets(line,512,stdin)) ) {

    if (n_files) {
       if (argv[i][0] == '-') continue;
       else strcpy(sac_file,argv[i]);
    } else {
       j=sscanf(line, "%s %f",sac_file,&t0);
       if (j==1) t0=0.;
    }
    fprintf(stderr,"%s\n",sac_file);

    if ( (trace=read_sac(sac_file,&hd)) == NULL 
	|| (nrdc == -2 && hd.a < -12340.)
	|| (rdc != 't' && hd.dist < -12340.)
        || (ntrace>0 && fabs(hd.delta-dt)>0.1*dt) ) {
	   fprintf(stderr,"error opening %s or missing/inconsistent head info\n",sac_file);
	   continue;
    }
    if (rdc == 't') 
       align = *((float *)(&hd) + 10 + nrdc) + t0 - hd.b;
    else
       align = hd.dist/reduce_vel - hd.b;

    ntrace++;
    if (ntrace == 1) {
       dt = hd.delta;
       ar = align;
       hd0=(SACHEAD *)malloc(sizeof(SACHEAD));
       assert(hd0 != NULL);
       memcpy(hd0, &hd, sizeof(SACHEAD));
       hd0->a = hd0->b + ar;
       if (wndw) {
          hd0->npts = rint((t2-t1)/dt);
	  hd0->b += ar+t1;
	  ar = -t1;
       }
       if (setT) {
          hd0->b -= hd0->a;
	  hd0->a = 0.;
       }
       if (unseteve) hd0->evla = hd0->evlo = -12345.;
       else hd0->stla = hd0->stlo = -12345.;
       hd0->dist = hd0->baz = hd0->user0 = 0.;
       nn = hd0->npts;
       src=(float *)malloc(nn*sizeof(float));
       assert(src != NULL);
       for(j=0;j<nn;j++) src[j] = 0.;
    }

    /* stacking */
    shift = rint((ar-align)/dt);
    start = shift;		if (start<0) start = 0;
    end   = hd.npts+shift;	if (end>nn) end = nn;
    if (wndw) {
       j = rint((align+t1)/dt); 	if (j>start) start=j;
       j = rint((align+t2)/dt)+shift;	if (j<end) end=j;
    }
    if (intg) cumsum(trace, hd.npts, dt);
    if (isqr) sqr(trace, hd.npts);
    if (smth) maver(trace, hd.npts, smth);
    maxA = 1.;
    if (norm) {
       maxA = 0.;
       for(j=start;j<end;j++) {
	  if (norm==2) maxA += trace[j-shift];
	  else if (trace[j-shift]>maxA) maxA = trace[j-shift];
       }
    }
    hd0->dist += hd.dist;
    hd0->baz += hd.baz;
    hd0->user0 += hd.user0;
    if (maxA>0.) {
	for (j=start; j<end; j++) src[j] += trace[j-shift]/maxA;
    }
    free(trace);

  }
  
  if (ntrace<1) return 0;

  ar = 1./ntrace;
  hd0->dist *= ar;
  hd0->gcarc = hd0->dist/111.2;
  hd0->baz *= ar;
  hd0->user0 *= ar;
  if (setbaz) {
     hd0->baz = baz;
     hd0->user0 = p;
  }
  for(j=0;j<nn;j++) src[j] *= ar;
  hd0->e = hd0->b + dt*hd0->npts;
  write_sac(outf, *hd0, src);
  
  return 0;

}
