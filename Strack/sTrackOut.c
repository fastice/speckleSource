#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "strack.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"
/*
   Output results
*/

static size_t fwriteOptionalBS1(void *ptr, size_t nitems, size_t size, FILE *fp, int32_t flags, int32_t byteOrder) {
   if(byteOrder == LSB) return fwrite(ptr, nitems, size, fp);
   else return fwriteBS(ptr, nitems, size, fp, flags);
}

// Obsolete ?
void sTrackOut(TrackParams *trackPar)
{
   FILE *fpR, *fpA, *fpC, *fpT;
   fprintf(stderr, "WRITING RESULTS\n\n");
   fpR = fopen(trackPar->outFileR, "w");
   fpA = fopen(trackPar->outFileA, "w");
   fpC = fopen(trackPar->outFileC, "w");
   fpT = fopen(trackPar->outFileT, "w");

   fwriteOptionalBS1(trackPar->offR[0], sizeof(float),
            (trackPar->nA) * trackPar->nR, fpR, FLOAT32FLAG, trackPar->byteOrder);
   fwriteOptionalBS1(trackPar->offA[0], sizeof(float),
            (trackPar->nA) * trackPar->nR, fpA, FLOAT32FLAG, trackPar->byteOrder);
   fwriteOptionalBS1(trackPar->corr[0], sizeof(float),
            (trackPar->nA) * trackPar->nR, fpC, FLOAT32FLAG, trackPar->byteOrder);
   fwriteOptionalBS1(trackPar->type[0], sizeof(char),
            (trackPar->nA) * trackPar->nR, fpT, BYTEFLAG, trackPar->byteOrder);

   fclose(fpR);
   fclose(fpA);
   fclose(fpC);
   fclose(fpT);
}

/*
  Write offsets, correlation, and match type to output file
 */
void writeOffsets(int32_t i, TrackParams *trackPar, FILE *fpR, FILE *fpA, FILE *fpC, FILE *fpT)
{
	fwriteOptionalBS1(trackPar->offR[i], sizeof(float), trackPar->nR, fpR, FLOAT32FLAG, trackPar->byteOrder);
	fflush(fpR);
	fwriteOptionalBS1(trackPar->offA[i], sizeof(float), trackPar->nR, fpA, FLOAT32FLAG, trackPar->byteOrder);
	fflush(fpA);
	fwriteOptionalBS1(trackPar->corr[i], sizeof(float), trackPar->nR, fpC, FLOAT32FLAG, trackPar->byteOrder);
	fflush(fpC);
	fwriteOptionalBS1(trackPar->type[i], sizeof(char), trackPar->nR, fpT, BYTEFLAG, trackPar->byteOrder);
	fflush(fpT);
}