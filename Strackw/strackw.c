#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "speckleSource/Strack/strack.h"
#include "strackw.h"
#include "math.h"

double SLat=-91.;
static void readArgs(int argc,char *argv[],char **parFile, int *floatFlag);
static void usage();
/* 
   Global variables definitions (NOT USED, NEED FOR LINKING ERS CODE
*/
int RangeSize=0;            /* Range size of complex image */
int AzimuthSize=0;          /* Azimuth size of complex image */
int BufferSize=0;           /* Size of nonoverlap region of the buffer */
int BufferLines =0;         /* # of lines of nonoverlap in buffer */
double RangePixelSize=0;    /* Range PixelSize */
double AzimuthPixelSize=0; /* Azimuth PixelSize */
int HemiSphere=0;
int DemType=0;
double Rotation=0;
char *Abuf1,*Abuf2,*Dbuf1,*Dbuf2;
int llConserveMem=0; /* Kluge to maintain backwards compat 9/13/06 */
float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
void *offBufSpace1,*offBufSpace2,*offBufSpace3,*offBufSpace4;
void *lBuf1,*lBuf2,*lBuf3,*lBuf4;

void main(int argc, char *argv[])
{   
	char *parFile;
	TrackParams trackPar;
	int noComplex;
	int floatFlag;
 	stateV sv1, sv2;
	
	readArgs(argc,argv,&parFile,&floatFlag);
	trackPar.floatFlag=floatFlag;
	trackPar.noComplex=TRUE;
	fprintf(stderr,"Using: %s\n", parFile);
	/*
	  Parse command file
	*/
	parseTrack(parFile, &trackPar);
	fprintf(stderr,"EdgePadR/A %i %i\n",trackPar.edgePadR,trackPar.edgePadA);
	/*
	  Parse image par files
	*/
	readOldPar(trackPar.parFile1,&(trackPar.imageP1),&sv1);
	readOldPar(trackPar.parFile2,&(trackPar.imageP2),&sv2);	
	/* 
	   Parse initial offsets
	*/
	parseInitialOffsets(&trackPar);
	/*
	  Get mask
	*/
	if(trackPar.maskFile != NULL) getMask(&trackPar);

	/*
	  Do matching
	*/
	corrTrackFast(&trackPar);
} 



static void readArgs(int argc,char *argv[],char **parFile,int *floatFlag)
{
	int n,i;
	char *argString;
	*floatFlag=TRUE;
	if( argc < 2 || argc > 3 ) usage();        /* Check number of args */ 
	*parFile = argv[argc-1];
	n = argc - 2;
	for(i=1; i <= n; i++) {
		argString = strchr(argv[i],'-');  
		if(strstr(argString,"integerComplex") != NULL)
			*floatFlag=FALSE;
		else usage();  
	}
	return;
}

 
static void usage()
{ 
	error("sTrack -integerComplex parFile \n");
}
