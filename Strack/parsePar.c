#include "clib/standard.h"
#include "stdio.h"
#include <stdlib.h>
#include <string.h>
#include "parfile.h"
#define MIN(a,b)  ( ( (a) < (b) ) ? (a) : (b) )

/*
    subroutine to read the SAR sensor parameters
    clw 24-dec-94
*/
  void parsePar(char *parFile, PROC_PAR *procp)
{
  int i;
  char *ret;
  int fret;
  FILE *procpf;

fprintf(stderr,"-- %s\n",parFile);
  procpf = openInputFile(parFile);
  rewind(procpf);			/* return to start of file */
  ret=fgets(procp->title,512,procpf);	/* title of the scene and processing     parameters */
  procp->title[strlen(procp->title)-1]='\0';/* strip off newline character */

  ret=fgets(procp->date,16,procpf);	/* date data acquired/generated (YYYY MM DD) */
  procp->date[strlen(procp->date)-1]='\0';	
  ret=fgets(procp->asciitime,24,procpf);	/* time of first record in signal data 
                                           file (HH MM SS.SSSS) */
  procp->asciitime[strlen(procp->asciitime)-1]='\0';	
  sscanf(procp->asciitime,"%d%d%lf",   &(procp->time.hour),&(procp->time.minute),&(procp->time.second));
fprintf(stderr,"%s",procp->asciitime);

  ret=fgets(procp->pol_ch,16,procpf);   /* polarization/channel to be processed out 
                                       of the set {HH, HV, VH, VV, CH1, CH2} */
  procp->pol_ch[strlen(procp->pol_ch)-1]='\0';

  fret=fscanf(procpf, "%lf %lf",&(procp->lat),&(procp->lon)); /* lat. and long. of
                                               scene center (decimal degrees) */

  fret=fscanf(procpf, "%lf",&(procp->track)); /* track angle of a/c (decimal deg) */

  fret=fscanf(procpf, "%lf",&(procp->alt));	/* altitude above geoid of the radar
                                           platform (m) */
  fret=fscanf(procpf, "%lf",&(procp->terra_alt)); /* avg height of terrain (m) */
    /* position of the SAR platform (m) at the center of  the presum data */	
  fret=fscanf(procpf, "%le %le %le",
	&(procp->pos.x),&(procp->pos.y),&(procp->pos.z));    
  fret=fscanf(procpf, "%lf %lf %lf ",&(procp->vel.x),&(procp->vel.y),
   &(procp->vel.z));  /* velocity vector (m/s), at center of the presum data */		
  fret=fscanf(procpf, "%lf %lf %lf ", &(procp->acc.x),&(procp->acc.y),
      &(procp->acc.z));  /* acc. vector (m/s^2), at center of the presum data */		
  fret=fscanf(procpf, "%lf",&(procp->prf));	/* radar prf (pulses/sec) */	
    /* doppler frequency polynomial coefficients fd=fdp[0]+fdp[1]*r+fdp[2]*r**2+
       fdp[3]*r**3 */
  fret=fscanf(procpf, "%le %le %le %le",&(procp->fdp[0]),&(procp->fdp[1]),
     &(procp->fdp[2]),&(procp->fdp[3])); 

  fret=fscanf(procpf, "%le",&(procp->td));	/* time delay between transmission of 
                                    pulse and first sample of echo (seconds) */
  fret=fscanf(procpf, "%d",&(procp->nr_offset)); /* offset (samples) to first range
           sample to process,each I/Q pair counts as 1 sample */	
  /* first,center,last slant range for the RAW SAR data */
  fret=fscanf(procpf, "%le %le %le",&(procp->r0),&(procp->r1),&(procp->r2));
  /* first,center,last slant range for the imaged swath */ 
  fret=fscanf(procpf, "%le %le %le",&(procp->ir0),&(procp->ir1),&(procp->ir2));
   /* slant range image pixel spacing, image range res.(m) */  
  fret=fscanf(procpf, "%lf %lf",&(procp->rpixsp),&(procp->ran_res)); 
  fret=fscanf(procpf, "%lf",&(procp->prfrac));/* fract of dop BW to process (Hz) */
  /* initial azimuth presumming factor (i.e. 4 for a 4:1 presum) */		
  fret=fscanf(procpf, "%d",&(procp->nprs_az));
 /* offset in echoes from the start of SAR signal data file */
  fret=fscanf(procpf, "%d",&(procp->loff));
  /* number of echoes to process/generate */	
  fret=fscanf(procpf, "%d",&(procp->nl));
  /* size of range FFT, # of range pixels to process, 2^N */
  fret=fscanf(procpf, "%d",&(procp->nrfft));	
  fret=fscanf(procpf, "%d",&(procp->nsub));		/* number of sub-images */
  /* #range looks, # azimuth looks/sub-image */
  fret=fscanf(procpf, "%d %d",&(procp->nlr),&(procp->nlaz));
  /* azimuth offset (m) of first image line from first prefiltered data line */ 
  fret=fscanf(procpf, "%lf",&(procp->azoff));
   /* azimuth image pixel spacing, azimuth image resolution */
  fret=fscanf(procpf, "%lf %lf",&(procp->azimsp), &(procp->azimres));
  /* image size in range, azimuth pixels */	
  fret=fscanf(procpf, "%d %d",&(procp->nrs),&(procp->nazs));	
  fret=fscanf(procpf, "%d", &(procp->nstate));		/* # state vectors */

  if (feof(procpf) != 0 || procp->nstate == 0){
    fprintf(stderr,
 "no orbit state vectors in the MSP processing parameter file: airborne SAR\n");

    procp->nstate = 0;				/* set # state vectors to 0 */ 
    rewind(procpf); 
    clearerr(procpf);
    return;
  }

  if(procp->nstate > MXST){
    fprintf(stderr,
    "\n# orbit state vectors in the MSP processing parameter file: %d\n",
    procp->nstate);
    fprintf(stderr,
    "WARNING: #state vectors exceeds structure allocation, increase MXST!\n");
  }


  /* time of first state vector, UTC seconds of day */
  fret=fscanf(procpf,"%lf",&(procp->t_state));
  /* time interval between state vectors (s) */  
  fret=fscanf(procpf,"%lf",&(procp->tis)); 

  for (i=0; i < MIN(procp->nstate,MXST); i++){
    fret=fscanf(procpf, "%le  %le  %le", 
      &(procp->state[i].pos.x), &(procp->state[i].pos.y), 
      &(procp->state[i].pos.z));  /* position */	 
    fret=fscanf(procpf, "%le  %le  %le", 
      &(procp->state[i].vel.x), &(procp->state[i].vel.y), 
      &(procp->state[i].vel.z));  /* velocity */
    procp->sTimes[i] = procp->t_state + i*procp->tis;
/*fprintf(stderr,"%f %f %f \n",procp->state[i].pos.x,procp->state[i].pos.y,procp->state[i].pos.z);*/
  } 
}

