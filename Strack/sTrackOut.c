#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "strack.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"
/* 
   Output results
*/

    void sTrackOut(TrackParams *trackPar)
{ 
   FILE *fpR, *fpA, *fpC, *fpT;  
  fprintf(stderr,"WRITING RESULTS\n\n");
  fpR = fopen(trackPar->outFileR,"w");
  fpA = fopen(trackPar->outFileA,"w");
  fpC = fopen(trackPar->outFileC,"w");
  fpT = fopen(trackPar->outFileT,"w");
  fwriteBS(trackPar->offR[0],sizeof(float),
     (trackPar->nA)*trackPar->nR,fpR,FLOAT32FLAG);
  fwriteBS(trackPar->offA[0],sizeof(float),
     (trackPar->nA)*trackPar->nR,fpA,FLOAT32FLAG);
  fwriteBS(trackPar->corr[0],sizeof(float),
     (trackPar->nA)*trackPar->nR,fpC,FLOAT32FLAG);
  fwriteBS(trackPar->type[0],sizeof(char),
     (trackPar->nA)*trackPar->nR,fpT,BYTEFLAG);

  fclose(fpR);
  fclose(fpA);
  fclose(fpC);
  fclose(fpT);
} 

