#include "stdio.h"
#include "stdint.h"
#include "string.h"
#include "clib/standard.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "cullst.h"
#include "math.h"
#include <stdlib.h>
#include <unistd.h>

#define LARGEINT -2e9
static char **mallocByteMat(int32_t nA, int32_t nR);
static void readArgs(int argc, char *argv[], CullParams *cullPar);
static void usage();
static void addSubtractSimOffsets(CullParams *cullPar, float mySign);
/*
   Global variables definitions
*/
double SLat = -91.;
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
	CullParams cullPar;
	FILE *fp;
	GDALAllRegister();

	readArgs(argc, argv, &cullPar);

	fprintf(stderr, "Infiles: \n %s\n %s\n %s\n %s\n %s\n", cullPar.inFileR, cullPar.inFileA, cullPar.inFileC, cullPar.inFileD, cullPar.inFileT);
	fprintf(stderr, "Outfiles: \n %s\n %s\n %s\n %s\n %s\n", cullPar.outFileR, cullPar.outFileA, cullPar.outFileC, cullPar.outFileD, cullPar.outFileT);
	fprintf(stderr, "cullIslandThresh %i\n", cullPar.islandThresh);
	/*
	   Load data
	*/
	loadCullData(&cullPar);
	/*
	   if useSim flag set, then read sim offsets, and subtract the intial guess prior to culling.
	 */
	if (cullPar.useSim == TRUE)
	{
		loadSimData(&cullPar);
		fprintf(stderr, "removing sim offsets \n");
		addSubtractSimOffsets(&cullPar, -1.0);
	}
	/*
	  load mask
	*/
	if (access("offsets.mask", F_OK) != -1 && cullPar.ignoreOffsets == FALSE)
	{
		/* read mask file if it is the same size, if not skip, and set flag FALSE */
		cullPar.maskFlag = loadCullMask(&cullPar);
	}
	else
		cullPar.maskFlag = FALSE;

	/*
	  cull data
	*/
	fprintf(stderr, "Cull data\n");
	cullSTData(&cullPar);
	cullSTData(&cullPar);
	cullSTData(&cullPar);
	/*
	  Compute stats for data
	*/
	fprintf(stderr, "Cull Stats\n");
	cullStats(&cullPar);
	/*
	   if useSim flag set, then add them back in after the residual difference was culled.
	 */
	if (cullPar.useSim == TRUE)
	{
		fprintf(stderr, "adding sim offsets back\n");
		addSubtractSimOffsets(&cullPar, 1.0);
	}
	/*
	  Compute stats for data
	*/
	fprintf(stderr, "Cull Smooth\n");
	cullSmooth(&cullPar);
	/*
	  cull islands
	*/
	if (cullPar.islandThresh > 0)
		cullIslands(&cullPar);

	/*
	  Output result
	*/
	fprintf(stderr, "Cull write\n");
	writeCullData(&cullPar);
}

static void addSubtractSimOffsets(CullParams *cullPar, float mySign)
{
	int32_t i, j;
	for (i = 0; i < cullPar->nR; i++)
		for (j = 0; j < cullPar->nA; j++)
		{
			if (cullPar->offA[j][i] < (-LARGEINT + 1) && cullPar->offSimA[j][i] < (-LARGEINT + 1))
			{
				cullPar->offA[j][i] += mySign * cullPar->offSimA[j][i];
				cullPar->offR[j][i] += mySign * cullPar->offSimR[j][i];
			}
		}
}

static void readArgs(int argc, char *argv[], CullParams *cullPar)
{
	char *argString;
	char *inBase, *outBase, *simBase;
	int32_t islandThresh;
	int32_t sr, sa;
	int32_t boxSize, nGood;
	float maxA, maxR, corrThresh;
	int32_t singleMT;
	int32_t i, n, sLen;

	if (argc < 3 || argc > 21)
	{
		fprintf(stderr, "to Many Argc %i %i\n", argc, 21);
		usage();
	} /* Check number of args */
	n = argc - 3;
	sr = 2;
	sa = 8;
	boxSize = 9;
	nGood = 17;
	maxA = 1.0;
	maxR = .75;
	singleMT = -1;
	islandThresh = -1;
	cullPar->ignoreOffsets = FALSE;
	cullPar->useSim = FALSE;
	corrThresh = 0.0;

	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "sr") != NULL)
		{
			sscanf(argv[i + 1], "%i", &sr);
			i++;
		}
		else if (strstr(argString, "sa") != NULL)
		{
			sscanf(argv[i + 1], "%i", &sa);
			i++;
		}
		else if (strstr(argString, "corrThresh") != NULL)
		{
			sscanf(argv[i + 1], "%f", &corrThresh);
			if(corrThresh < 0 || corrThresh >=1) 
			{
				fprintf(stderr, "\n\ncoorThresh must be in interval (0,1)\n\n");
				usage();
			}
			i++;
		}
		else if (strstr(argString, "boxSize") != NULL)
		{
			sscanf(argv[i + 1], "%i", &boxSize);
			i++;
		}
		else if (strstr(argString, "nGood") != NULL)
		{
			sscanf(argv[i + 1], "%i", &nGood);
			i++;
		}
		else if (strstr(argString, "singleMT") != NULL)
		{
			sscanf(argv[i + 1], "%i", &singleMT);
			if (singleMT < 1 || singleMT > 3)
			{
				fprintf(stderr, "singleMT should be in range 0 to 3\n");
				usage();
			}
			i++;
		}
		else if (strstr(argString, "maxA") != NULL)
		{
			sscanf(argv[i + 1], "%f", &maxA);
			i++;
		}
		else if (strstr(argString, "maxR") != NULL)
		{
			sscanf(argv[i + 1], "%f", &maxR);
			i++;
		}
		else if (strstr(argString, "islandThresh") != NULL)
		{
			sscanf(argv[i + 1], "%i", &islandThresh);
			i++;
		}
		else if (strstr(argString, "useSim") != NULL)
		{
			cullPar->useSim = TRUE;
			simBase = argv[i + 1];
			i++;
		}
		else if (strstr(argString, "ignoreOffsets") != NULL)
		{
			cullPar->ignoreOffsets = TRUE;
		}
		else
		{
			fprintf(stderr, "%i %s\n", i, argv[i]);
			usage();
		}
	}
	cullPar->islandThresh = islandThresh;
	cullPar->sR = sr;
	cullPar->sA = sa;
	cullPar->bR = boxSize;
	cullPar->bA = boxSize;
	cullPar->nGood = nGood;
	cullPar->maxA = maxA;
	cullPar->maxR = maxR;
	cullPar->corrThresh = corrThresh;
	cullPar->singleMT = singleMT;

	inBase = argv[argc - 2];
	outBase = argv[argc - 1];

	fprintf(stderr, "|%s|\n", inBase);
	fprintf(stderr, "|%s|\n", outBase);

	const char *baseNames[] = {outBase, outBase, inBase, inBase, outBase, outBase, outBase, outBase, outBase};
	const char *exts[] = {".dr", ".da",".sr", ".sa", ".cc", ".mt", ".dat", ".vrt", ".mt.vrt"};
	const int numExts = sizeof(exts) / sizeof(exts[0]);
	sLen = strlen(outBase);

	for (int i = 0; i < numExts; i++)
	{
		const size_t len = sLen + strlen(exts[i]) + 1;
		char *filename = (char *)malloc(sizeof(char) * len);
		filename[0] = '\0';
		strcpy(filename, baseNames[i]);
		strcat(filename, exts[i]);

		switch (i)
		{
		case 0:
			cullPar->outFileR = filename;
			break;
		case 1:
			cullPar->outFileA = filename;
			break;
		case 2:
			cullPar->outFileSR = filename;
			break;
		case 3:
			cullPar->outFileSA = filename;
			break;
		case 4:
			cullPar->outFileC = filename;
			break;
		case 5:
			cullPar->outFileT = filename;
			break;
		case 6:
			cullPar->outFileD = filename;
			break;
		case 7:
			cullPar->outFileVRT = filename;
			break;
		case 8:
			cullPar->outFileVRTMT = filename;
			break;
		}
	}
	const char *extsI[] = {".dr", ".da", ".cc", ".mt", ".dat", ".vrt", ".mt.vrt"};
	const int numExtsI = sizeof(extsI) / sizeof(extsI[0]);
	sLen = strlen(inBase);
	for (int i = 0; i < numExtsI; i++)
	{
		const size_t len = sLen + strlen(extsI[i]) + 1;
		char *filename = (char *)malloc(sizeof(char) * len);
		filename[0] = '\0';
		strcpy(filename, inBase);
		strcat(filename, extsI[i]);

		switch (i)
		{
		case 0:
			cullPar->inFileR = filename;
			break;
		case 1:
			cullPar->inFileA = filename;
			break;
		case 2:
			cullPar->inFileC = filename;
			break;
		case 3:
			cullPar->inFileT = filename;
			break;
		case 4:
			cullPar->inFileD = filename;
			break;
		case 5:
			cullPar->inFileVRT = filename;
			break;
		case 6:
			cullPar->inFileVRTMT = filename;
			break;
		}
	}

	if (cullPar->useSim == TRUE)
	{
		const char *exts[] = {".dr", ".da", ".dat", ".vrt", ".cc"};
		const int numExts = sizeof(exts) / sizeof(exts[0]);
		sLen = strlen(simBase);
		for (int i = 0; i < numExts; i++)
		{
			const size_t len = sLen + strlen(exts[i]) + 1;
			char *filename = (char *)malloc(sizeof(char) * len);
			filename[0] = '\0';
			strcpy(filename, simBase);
			strcat(filename, exts[i]);

			switch (i)
			{
			case 0:
				cullPar->inFileSimR = filename;
				break;
			case 1:
				cullPar->inFileSimA = filename;
				break;
			case 2:
				cullPar->inFileSimD = filename;
				break;
			case 3:
				cullPar->inFileSimVRT = filename;
				break;
			}
		}
	}
	cullPar->metaData = NULL;
	cullPar->metaDataMT = NULL;

	return;
}

static void usage()
{
	error("cullst -useSim offsets -ignoreOffsets -islandThresh islandThresh -singleMT singleMT -maxR maxR -maxA maxA -nGood nGood -boxSize boxSize -sr sr -sa sa inbase outbase\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n%s\n",
		  "where",
		  "useSim is a simulated offsets file to subtract off in cull",
		  "sr,sa = smoothing window width (full) to apply (if odd then 3->111, if even then 4-> 0.5 1 1 1 0.5)",
		  "islandThresh = cull isolated islands < islandThresh in diameter",
		  "singleMT = only retain this match type (1,2,or3)",
		  "maxR,maxA = max deviation from local median",
		  "corrThresh = min correlation (default 0 to pass all)",
		  "ignoreOffsets = do not use information from offset file to cull large matches");
}
