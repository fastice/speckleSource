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
char *Abuf1, *Abuf2, *Dbuf1, *Dbuf2;
int32_t llConserveMem = 0;			/* Klugxse to maintain backwards compat 9/13/06 */
float *AImageBuffer, *DImageBuffer; /* Kluge 05/31/07 to seperate image buffers */
void *offBufSpace1, *offBufSpace2, *offBufSpace3, *offBufSpace4;
void *lBuf1, *lBuf2, *lBuf3, *lBuf4;

int main(int argc, char *argv[])
{
	char *parFile;
	TrackParams trackPar;
	int32_t noComplex;
	int32_t floatFlag;
	stateV sv1, sv2;
	int32_t byteOrder;
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
	/*
	  Parse image par files
	*/
	readOldPar(trackPar.parFile1, &(trackPar.imageP1), &sv1);
	readOldPar(trackPar.parFile2, &(trackPar.imageP2), &sv2);
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
