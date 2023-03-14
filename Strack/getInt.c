#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "strack.h"
#include "math.h"
#include "sys/stat.h"
/*
  Malloc space, read params for, and input complex interferogram
*/

void getInt(TrackParams *trackPar)
{   
	inputImageStructure inputImage;
	fftw_complex *buf;
	size_t bufSize;
	strackComplex *inBuf;
	FILE *fp;
	struct stat st;
	int i,j, nLines;

	fprintf(stderr,"get int \n");

	if(trackPar->intGeodat==NULL) error("getInt: Missing geodat filename");

	inputImage.stateFlag=TRUE;
	parseInputFile(trackPar->intGeodat, &inputImage);
	trackPar->latc=inputImage.latControlPoints[0];
	trackPar->intDat.nr = inputImage.rangeSize; 
	trackPar->intDat.na = inputImage.azimuthSize; 
	trackPar->intDat.nrl = inputImage.nRangeLooks; 
	trackPar->intDat.nal= inputImage.nAzimuthLooks; 
	trackPar->lambda=inputImage.par.lambda;
	inBuf = (strackComplex *) malloc(trackPar->intDat.nr *sizeof(strackComplex));       
	trackPar->intDat.intf = (fftw_complex **) malloc(trackPar->intDat.na * sizeof(fftw_complex *));
	bufSize =(long)trackPar->intDat.na * (long)trackPar->intDat.nr *(long)sizeof(fftw_complex); 
	buf = (fftw_complex *) malloc(bufSize);
	if(buf == NULL) {
		fprintf(stderr,"Cannot malloc interferogram buffer, proceeding with no interferogram\n");
		trackPar->intDat.intf = NULL;
		return;
	}
	fprintf(stderr,"Reading Interferogram nr,na, bufsize %d %d %u \n", trackPar->intDat.nr,  trackPar->intDat.na, bufSize);

	for(i=0; i < trackPar->intDat.na; i++) {
		trackPar->intDat.intf[i] = &(buf[i*trackPar->intDat.nr]);
		/* Fill buffer with something close to zero */
		for(j=0; j < trackPar->intDat.nr; j++) {
			trackPar->intDat.intf[i][j].re = 1;
			trackPar->intDat.intf[i][j].im = 1;			
		}
	}
       
	/* Allow interferogram to be smaller than image - only read in part with data */
	stat(trackPar->intFile,  &st);
	/* Actual number of valid lines */
	nLines = st.st_size/(trackPar->intDat.nr * 4);
	fprintf(stderr,"geodat lines %i int %i\n",trackPar->intDat.na, nLines);
	fp = fopen(trackPar->intFile,"r");
	/*
	  Read interferogram - if short the last part will have fill values defined above.
	 */
	for(i=0; i < nLines; i++) {
		freadBS(inBuf,sizeof(strackComplex),trackPar->intDat.nr,fp,INT16FLAG);
		for(j=0; j < trackPar->intDat.nr; j++) {
			trackPar->intDat.intf[i][j].re = (float)inBuf[j].r;
			trackPar->intDat.intf[i][j].im = (float)inBuf[j].i;
		}
	}
	fclose(fp);
	free(inBuf);
} 

