#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "speckleSource/Strack/strack.h"
#include "strackw.h"
#include "math.h"

double SLat = -91.;
static void readArgs(int argc, char *argv[], char **parFile, int32_t *floatFlag, int32_t *byteOrder);
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

int32_t llConserveMem = 0;			/* Klugxse to maintain backwards compat 9/13/06 */


int main(int argc, char *argv[])
{
	char *parFile;
	TrackParams trackPar;
	int32_t noComplex;
	int32_t floatFlag;
	stateV sv1, sv2;
	int32_t byteOrder;
	GDALDatasetH hDS1, hDS2;
	GDALAllRegister();
	readArgs(argc, argv, &parFile, &floatFlag, &byteOrder);
	trackPar.floatFlag = floatFlag;
	trackPar.noComplex = TRUE;
	trackPar.byteOrder = byteOrder;
	fprintf(stderr, "Using: %s\n", parFile);
	/*
	  Parse command file
	*/
	parseTrack(parFile, &trackPar);
	fprintf(stderr, "EdgePadR/A %i %i\n", trackPar.edgePadR, trackPar.edgePadA);
	// Parse First Image
	fprintf(stderr, "%s  %s\n",trackPar.vrtFile1, trackPar.vrtFile2 );

	if(trackPar.parFile1 != NULL) 
	{
		readOldPar(trackPar.parFile1, &(trackPar.imageP1), &sv1);
		trackPar.hBand1 = NULL;
	} 
	else if(trackPar.vrtFile1 != NULL)
	{
		fprintf(stderr, "Parsing VRT file for image 1\n");
		parseSLCVrtNew(trackPar.vrtFile1, &(trackPar.imageP1), &sv1, &byteOrder, &hDS1, &trackPar.hBand1);
		trackPar.fpI1 = NULL;
		//trackPar.hBand1 = NULL;
	} 
	else error("Could not read par or vrt file for image 1\n");

	// Parse Second Image
	if(trackPar.parFile2 != NULL) 
	{
		readOldPar(trackPar.parFile2, &(trackPar.imageP2), &sv2);
		trackPar.hBand2 = NULL;
	} 
	else if(trackPar.vrtFile2 != NULL)
	{
		parseSLCVrtNew(trackPar.vrtFile2, &(trackPar.imageP2), &sv2, &byteOrder, &hDS2, &trackPar.hBand2);
		trackPar.fpI2 = NULL;
	}
	else error("Could not read par or vrt file for image 1\n");


	if(trackPar.parFile1 != NULL) 
	{
		readOldPar(trackPar.parFile1, &(trackPar.imageP1), &sv1);
		trackPar.hBand1 = NULL;
	} 
	else if(trackPar.vrtFile1 != NULL)
	{
		fprintf(stderr, "Parsing VRT file for image 1\n");
		parseSLCVrtNew(trackPar.vrtFile1, &(trackPar.imageP1), &sv1, &byteOrder, &hDS1, &trackPar.hBand1);
		trackPar.fpI1 = NULL;
		fprintf(stderr, "finished parsing VRT file for image 1\n");
		//parseSLCVrt(trackPar.vrtFile1, &(trackPar.imageP1), &sv1, &byteOrder);
	} 
	else error("Could not read par or vrt file for image 1\n");

	// Parse Second Image
	if(trackPar.parFile2 != NULL) 
	{
		readOldPar(trackPar.parFile2, &(trackPar.imageP2), &sv2);
		trackPar.hBand2 = NULL;
	} 
	else if(trackPar.vrtFile2 != NULL)
	{
		fprintf(stderr, "Parsing VRT file for image 2\n");
		parseSLCVrtNew(trackPar.vrtFile2, &(trackPar.imageP2), &sv2, &byteOrder, &hDS2, &trackPar.hBand2);
		trackPar.fpI2 = NULL;
		fprintf(stderr, "finished parsing VRT file for image 2\n");
		//parseSLCVrt(trackPar.vrtFile2, &(trackPar.imageP2), &sv2, &byteOrder);
	}
	else error("Could not read par or vrt file for image 1\n");
	// Override command line/default if vrt specified 
	if(byteOrder >= 0) trackPar.byteOrder = byteOrder;
	//
	fprintf(stderr, "EdgePadR/A %i %i\n", trackPar.edgePadR, trackPar.edgePadA);
	/*
	   Parse initial offsets
	*/
	parseInitialOffsets(&trackPar);
	fprintf(stderr, "EdgePadR/A %i %i\n", trackPar.edgePadR, trackPar.edgePadA);
	/*
	  Get mask
	*/
	if (trackPar.maskFile != NULL)
		getMask(&trackPar);
	/*
	  Do matching
	*/
fprintf(stderr, "EdgePadR/A %i %i\n", trackPar.edgePadR, trackPar.edgePadA);
	corrTrackFast(&trackPar);
}

static void readArgs(int argc, char *argv[], char **parFile, int32_t *floatFlag, int32_t *byteOrder)
{
	int32_t n, i;
	char *argString;
	*floatFlag = TRUE;
	*byteOrder = MSB;
	if (argc < 2 || argc > 4)
		usage(); /* Check number of args */
	*parFile = argv[argc - 1];
	n = argc - 2;
	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "integerComplex") != NULL)
			*floatFlag = FALSE;
		else if (strstr(argString, "-LSB") != NULL)
			*byteOrder = LSB;
		else
			usage();
	}
	return;
}

static void usage()
{
	error("sTrack -integerComplex -LSB parFile \n\tLSB  use LSB for both output and flaoting point input\n");
}
