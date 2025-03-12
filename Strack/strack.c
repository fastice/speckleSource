#include "string.h"
#include "strack.h"
#include "math.h"

double SLat = -91.;

static void readArgs(int32_t argc, char *argv[], char **parFile, int32_t *noComplex, int32_t *floatFlag, 
	int32_t *hanningFlag, int32_t *legacyFlag, int32_t *gaussFlag, int32_t *maxTries, int32_t *byteOrder);
static void usage();
/*
   Global variables definitions (NOT USED, NEED FOR LINKING ERS CODE
*/
int32_t RangeSize = 0;		 /* Range size of complex image */
int32_t AzimuthSize = 0;	 /* Azimuth size of complex image */
int32_t BufferSize = 0;		 /* Size of nonoverlap region of the buffer */
int32_t BufferLines = 0;	 /* # of lines of nonoverlap in buffer */
double RangePixelSize = 0;	 /* Range PixelSize */
double AzimuthPixelSize = 0; /* Azimuth PixelSize */
int32_t HemiSphere = 0;
int32_t DemType = 0;
double Rotation = 0;
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
int32_t llConserveMem = 0;			/* Kluge to maintain backwards compat 9/13/06 */
float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
void *lBuf1, *lBuf2, *lBuf3, *lBuf4;

int main(int argc, char *argv[])
{
	char *parFile;
	TrackParams trackPar;
	int32_t noComplex;
	int32_t floatFlag, maxTries;
	inputImageStructure inputImage;
	int32_t hanningFlag, legacyFlag, gaussFlag, byteOrder;
	stateV sv1, sv2;
	readArgs(argc, argv, &parFile, &noComplex, &floatFlag, &hanningFlag, &legacyFlag, &gaussFlag, &maxTries, &byteOrder);
	trackPar.floatFlag = floatFlag;
	trackPar.noComplex = noComplex;
	trackPar.hanningFlag = hanningFlag;
	trackPar.gaussFlag = gaussFlag;
	trackPar.legacyFlag = legacyFlag; /* Flag to use stuff for RADARSAT - doppler etc */
	trackPar.maxTries = maxTries;
	trackPar.byteOrder = byteOrder;
	fprintf(stderr,"Byte order %i  [LSB %i, MSB %i]\n", trackPar.byteOrder, LSB, MSB);
	fprintf(stderr, "maxTries = %i\n", maxTries);
	GDALAllRegister();
	if (noComplex == TRUE)
	{
		fprintf(stderr, "\n\n*** noComplex flag set- amplitude matching only ***\n\n");
	}
	fprintf(stderr, "Using: %s\n", parFile);
	/*
	  Parse command file
	*/
	parseTrack(parFile, &trackPar);

	/*
	  Parse image par files
	*/
	readOldPar(trackPar.parFile1, &(trackPar.imageP1), &sv1);
	readOldPar(trackPar.parFile2, &(trackPar.imageP2), &sv2);
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
	if (trackPar.intFile != NULL)
	{
		if (trackPar.noComplex == FALSE)
		{
			getInt(&trackPar);
			fprintf(stderr, "Lambda %f\n", trackPar.lambda);
		}
		if (trackPar.intDat.intf == NULL)
		{ /* If not enough mem for int, don't use, and proceed without */
			trackPar.intFile = NULL;
			trackPar.intFlag = FALSE;
		}
	}
	/* Case with no int32_t by design, or for lack of memory */
	if (trackPar.intFile == NULL)
	{ /* Added 8/8/2016 for baseline only correction */
		if (trackPar.intGeodat == NULL)
			error("getInt: Missing geodat filename");
		/*if(trackPar.noComplex==FALSE) {*/
		inputImage.stateFlag = TRUE;
		parseInputFile(trackPar.intGeodat, &inputImage);
		trackPar.lambda = inputImage.par.lambda;
		trackPar.latc = inputImage.latControlPoints[0];
		/*}*/
	}
	/*
	  Get mask
	*/
	if (trackPar.maskFile != NULL)
		getMask(&trackPar);
	fprintf(stderr, "C %s\n", trackPar.maskFile);
	/*
	  Do matching
	*/
	speckleTrack(&trackPar);
}

static void readArgs(int32_t argc, char *argv[], char **parFile, int32_t *noComplex,
					 int32_t *floatFlag, int32_t *hanningFlag, int32_t *legacyFlag, 
					 int32_t *gaussFlag, int32_t *maxTries, int *byteOrder)
{
	int32_t n, i;
	char *argString;
	*floatFlag = TRUE;
	if (argc < 2 || argc > 5)
		usage(); /* Check number of args */
	*parFile = argv[argc - 1];
	n = argc - 2;
	*noComplex = FALSE;
	*hanningFlag = TRUE;
	*legacyFlag = FALSE;
	*gaussFlag = FALSE;
	*byteOrder = MSB;
	*maxTries = 2;
	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "noComplex") != NULL)
			*noComplex = TRUE;
		else if (strstr(argString, "integerComplex") != NULL)
			*floatFlag = FALSE;
		else if (strstr(argString, "-noHanning") != NULL)
			*hanningFlag = FALSE;
		else if (strstr(argString, "-legacy") != NULL)
			*legacyFlag = TRUE;
		else if (strstr(argString, "-gauss") != NULL)
			*gaussFlag = TRUE;
		else if (strstr(argString, "-singleAmp") != NULL)
			*maxTries = 1;
		else if (strstr(argString, "-LSB") != NULL)
			*byteOrder = LSB;
		else
			usage();
	}
	return;
}

static void usage()
{
	error("sTrack -noComplex -legacy -singleAmp -gauss -integerComplex -noHanning -LSB parFile \n\tLSB  use LSB for both output and flaoting point input\n ");
}
