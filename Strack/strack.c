#include"string.h"
#include "strack.h"
#include "math.h"


static void readArgs(int argc,char *argv[],char **parFile,int *noComplex, int *floatFlag, int *hanningFlag);
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
	inputImageStructure inputImage;
	int hanningFlag;
	stateV sv1, sv2;
	readArgs(argc,argv,&parFile,&noComplex,&floatFlag, &hanningFlag);
	trackPar.floatFlag=floatFlag;
	trackPar.noComplex=noComplex;
	trackPar.hanningFlag=hanningFlag;
	if(noComplex==TRUE) {
		fprintf(stderr, "\n\n*** noComplex flag set- amplitude matching only ***\n\n");
	}
	fprintf(stderr,"Using: %s\n", parFile);
	/*
	  Parse command file
	*/
	parseTrack(parFile, &trackPar);

	/*
	  Parse image par files
	*/
	readOldPar(trackPar.parFile1,&(trackPar.imageP1),&sv1);
	readOldPar(trackPar.parFile2,&(trackPar.imageP2),&sv2);
	/* 
	   Parse baseline
	*/
	parseBase(&trackPar);
	/* 
	   Parse initial offsets
	   if(trackPar.polyShift==TRUE)
	   removed if 6/2/04 to use offset polynomial as backup
	*/
	parseInitialOffsets(&trackPar);
	/*
	  Get inteferogram
	*/

	if(trackPar.intFile != NULL) {
		if(trackPar.noComplex == FALSE) getInt(&trackPar);
		fprintf(stderr,"Lambda %f\n",trackPar.lambda);
	} else { /* Added 8/8/2016 for baseline only correction */
		if(trackPar.intGeodat==NULL ) error("getInt: Missing geodat filename");
		/*if(trackPar.noComplex==FALSE) {*/
			inputImage.stateFlag=TRUE;
			parseInputFile(trackPar.intGeodat, &inputImage);
			trackPar.lambda=inputImage.par.lambda;
			trackPar.latc=inputImage.latControlPoints[0];
		/*}*/
	}
	/*
	  Get mask
	*/
	if(trackPar.maskFile != NULL) getMask(&trackPar);
	fprintf(stderr,"C %s\n",trackPar.maskFile);
	/*
	  Do matching
	*/
	speckleTrack(&trackPar);
} 



static void readArgs(int argc,char *argv[],char **parFile,int *noComplex,int *floatFlag, int *hanningFlag)
{
	int n,i;
	char *argString;
	*floatFlag=TRUE;
	if( argc < 2 || argc > 4 ) usage();        /* Check number of args */ 
	*parFile = argv[argc-1];
	n = argc - 2;
	*noComplex=FALSE;
	*hanningFlag=TRUE;
	for(i=1; i <= n; i++) {
		argString = strchr(argv[i],'-');  
		if(strstr(argString,"noComplex") != NULL)
			*noComplex=TRUE;
		else if(strstr(argString,"integerComplex") != NULL)
			*floatFlag=FALSE;
		else if(strstr(argString,"-noHanning") != NULL)
		       *hanningFlag=FALSE;
		else usage();  
	}
	return;
}

 
static void usage()
{ 
	error("sTrack -noComplex -integerComplex -noHanning parFile \n");
}
