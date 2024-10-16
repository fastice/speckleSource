#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "strack.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"

static void parseValues(RDF *commands, char *keyword, char *format,
						void *d1, void *d2, void *d3, int32_t n);

void printTrackPar(TrackParams *trackPar);

void parseTrack(char *parFile, TrackParams *trackPar)
{
	RDF *rdfParams;
	void *d2, *d3;
	char *tmp;
	int32_t sLen;
	char buf[2048];
	d2 = NULL;
	d3 = NULL;
	rdfParams = NULL;
	rdfParams = rdfParse(parFile, rdfParams);
	/*
	  Image files
	*/
	if ((trackPar->imageFile1 = rdfValue(rdfParams, "image1")) == NULL)
		error("parseTrack: Missing image1");

	if ((trackPar->imageFile2 = rdfValue(rdfParams, "image2")) == NULL)
		error("parseTrack: Missing image2");

	if ((trackPar->parFile1 = rdfValue(rdfParams, "image1par")) == NULL)
		error("parseTrack: Missing image1par");

	if ((trackPar->parFile2 = rdfValue(rdfParams, "image2par")) == NULL)
		error("parseTrack: Missing image2par");
	/* Read parameters for mask file */
	trackPar->maskFlag = FALSE;
	if ((trackPar->maskFile = rdfValue(rdfParams, "maskfile")) != NULL)
	{
		trackPar->maskFlag = TRUE;
		trackPar->maskType = GEODATMASK;
		if ((trackPar->maskGeodat = rdfValue(rdfParams, "maskgeodat")) == NULL)
			error("parseTrack: Missing maskgeodat");
		fprintf(stderr, "\n***********************************************\n");
		fprintf(stderr, "Using maskfile %s\n", trackPar->maskFile);
		fprintf(stderr, "\n***********************************************\n");
	}
	else
	{
		if ((trackPar->maskFile = rdfValue(rdfParams, "offsetmaskfile")) != NULL)
		{
			trackPar->maskFlag = TRUE;
			trackPar->maskType = OFFSETMASK;
			if ((trackPar->maskGeodat = rdfValue(rdfParams, "offsetmaskdat")) == NULL)
			{
				buf[0] = '\0';
				strcpy(buf, trackPar->maskFile);
				strcat(buf, ".dat");
			}
			else
				fprintf(stderr, "\n***********************************************\n");
			fprintf(stderr, "Using maskfile %s\n", trackPar->maskFile);
			fprintf(stderr, "Using maskfile data %s\n", trackPar->maskGeodat);
			fprintf(stderr, "\n***********************************************\n");
		}
		trackPar->maskvrt = rdfValue(rdfParams, "offsetmaskvrt");
	}

	if ((trackPar->initialOffsetFile = rdfValue(rdfParams, "initialoffsetfile")) == NULL)
		error("parseTrack: Missing initialoffsetfile");

	trackPar->intGeodat = NULL;
	trackPar->baseParams = NULL;
	if ((trackPar->intFile = rdfValue(rdfParams, "intfile")) == NULL)
		trackPar->intFlag = FALSE;
	else
	{
		trackPar->intFlag = TRUE;
		if ((trackPar->intGeodat = rdfValue(rdfParams, "intgeodat")) == NULL)
			error("parseTrack: Missing intgeodat");
	}
	if ((trackPar->baseParams = rdfValue(rdfParams, "baseparams")) == NULL)
		trackPar->baseFlag = FALSE;
	else
	{
		trackPar->baseFlag = TRUE;
		fprintf(stderr, "Using Baseline %s\n", trackPar->baseParams);
		if ((trackPar->intGeodat = rdfValue(rdfParams, "intgeodat")) == NULL)
			error("parseTrack: Missing intgeodat");
	}
	if (trackPar->noComplex == TRUE)
		trackPar->baseFlag = FALSE;
	if ((trackPar->initOffsetsFile = rdfValue(rdfParams, "initShift")) == NULL)
	{
		trackPar->polyShift = TRUE;
		fprintf(stderr, "Polyshift = TRUE\n");
	}
	else
	{
		trackPar->polyShift = FALSE;
		fprintf(stderr, "Polyshift = FALSE\n");
		fprintf(stderr, "Using initial shifts from %s\n", trackPar->initOffsetsFile);
	}
	tmp = rdfValue(rdfParams, "outputfile");
	if (tmp != NULL)
	{
		sLen = strlen(tmp);
		trackPar->outFileR = (char *)malloc(sizeof(char) * sLen + 4);
		trackPar->outFileA = (char *)malloc(sizeof(char) * sLen + 4);
		trackPar->outFileC = (char *)malloc(sizeof(char) * sLen + 4);
		trackPar->outFileT = (char *)malloc(sizeof(char) * sLen + 4);
		trackPar->outFileD = (char *)malloc(sizeof(char) * sLen + 5);
		trackPar->vrtFile = (char *)malloc(sizeof(char) * sLen + 5);
		trackPar->MTvrtFile = (char *)malloc(sizeof(char) * sLen + 8);
		strcpy(trackPar->outFileR, tmp);
		strcpy(trackPar->outFileA, tmp);
		strcpy(trackPar->outFileC, tmp);
		strcpy(trackPar->outFileT, tmp);
		strcpy(trackPar->outFileD, tmp);
		strcpy(trackPar->vrtFile, tmp);
		strcpy(trackPar->MTvrtFile, tmp);
		strcat(trackPar->outFileR, ".dr");
		strcat(trackPar->outFileA, ".da");
		strcat(trackPar->outFileC, ".cc");
		strcat(trackPar->outFileT, ".mt");
		strcat(trackPar->outFileD, ".dat");
		strcat(trackPar->vrtFile, ".vrt");
		strcat(trackPar->MTvrtFile, ".mt.vrt");
		fprintf(stderr, "Outfiles: %s\n%s\n%s\n%s\n%s\n%s\n",
				trackPar->outFileR, trackPar->outFileA, trackPar->outFileC,
				trackPar->outFileR, trackPar->outFileA, trackPar->vrtFile);
	}
	else
		error("parseTrack: missing ouput file");

	parseValues(rdfParams, "rstart", "%i", &(trackPar->rStart), d2, d3, 1);
	parseValues(rdfParams, "astart", "%i", &(trackPar->aStart), d2, d3, 1);
	parseValues(rdfParams, "deltar", "%i", &(trackPar->deltaR), d2, d3, 1);
	parseValues(rdfParams, "deltaa", "%i", &(trackPar->deltaA), d2, d3, 1);

	parseValues(rdfParams, "nR", "%i", &(trackPar->nR), d2, d3, 1);
	parseValues(rdfParams, "nA", "%i", &(trackPar->nA), d2, d3, 1);

	parseValues(rdfParams, "wr", "%i", &(trackPar->wR), d2, d3, 1);
	parseValues(rdfParams, "wa", "%i", &(trackPar->wA), d2, d3, 1);

	parseValues(rdfParams, "wra", "%i", &(trackPar->wRa), d2, d3, 1);
	parseValues(rdfParams, "waa", "%i", &(trackPar->wAa), d2, d3, 1);

	trackPar->edgePad = 0;
	if (rdfValue(rdfParams, "edgepad") != NULL)
	{
		if (sscanf(rdfValue(rdfParams, "edgepad"), "%i", &(trackPar->edgePad)) != 1)
			trackPar->edgePad = 0;
	}
	trackPar->edgePadR = trackPar->edgePad;
	trackPar->edgePadA = trackPar->edgePad;

	trackPar->scaleFactor = 1;
	if (rdfValue(rdfParams, "scalefactor") != NULL)
	{
		if (sscanf(rdfValue(rdfParams, "scalefactor"), "%i", &(trackPar->scaleFactor)) != 1)
			trackPar->scaleFactor = 1;
	}
	fprintf(stderr, "Scale Factor %i\n", trackPar->scaleFactor);

	if (rdfValue(rdfParams, "edgeR") != NULL)
	{
		if (sscanf(rdfValue(rdfParams, "edgeR"), "%i", &(trackPar->edgePadR)) != 1)
			trackPar->edgePadR = trackPar->edgePad;
	}

	if (rdfValue(rdfParams, "edgeA") != NULL)
	{
		if (sscanf(rdfValue(rdfParams, "edgeA"), "%i", &(trackPar->edgePadA)) != 1)
			trackPar->edgePadR = trackPar->edgePad;
	}

	fprintf(stderr, "edgepad = %i %i %i\n", trackPar->edgePad, trackPar->edgePadR, trackPar->edgePadA);

	parseValues(rdfParams, "navgR", "%i", &(trackPar->navgR), d2, d3, 1);
	parseValues(rdfParams, "navgA", "%i", &(trackPar->navgA), d2, d3, 1);

	trackPar->rStart /= trackPar->scaleFactor;
	trackPar->aStart /= trackPar->scaleFactor;
	trackPar->deltaR /= trackPar->scaleFactor;
	trackPar->deltaA /= trackPar->scaleFactor;
	printTrackPar(trackPar);
}

/*
  Debug routine
*/
void printTrackPar(TrackParams *trackPar)
{
	fprintf(stderr, "\nimageFile1 = %s\n", trackPar->imageFile1);
	fprintf(stderr, "imageFile2 = %s\n", trackPar->imageFile2);
	fprintf(stderr, "parFile1 = %s\n", trackPar->parFile1);
	fprintf(stderr, "parFile2 = %s\n", trackPar->parFile2);
	fprintf(stderr, "initialOffsetFile = %s\n", trackPar->initialOffsetFile);
	fprintf(stderr, "1\n");
	if (trackPar->intFile != NULL)
		fprintf(stderr, "intFile = %s\n", trackPar->intFile);
	fprintf(stderr, "2\n");
	if (trackPar->intGeodat != NULL)
		fprintf(stderr, "intgeodat = %s\n", trackPar->intGeodat);
	fprintf(stderr, "3\n");
	if (trackPar->baseParams != NULL)
		fprintf(stderr, "baseParams = %s\n", trackPar->baseParams);

	fprintf(stderr, "\nscaleFactor = %i\n", trackPar->scaleFactor);
	fprintf(stderr, "\nrstart = %i\n", trackPar->rStart);
	fprintf(stderr, "astart = %i\n", trackPar->aStart);
	fprintf(stderr, "deltar = %i\n", trackPar->deltaR);
	fprintf(stderr, "deltaa = %i\n", trackPar->deltaA);

	fprintf(stderr, "\nnR = %i\n", trackPar->nR);
	fprintf(stderr, "nA = %i\n", trackPar->nA);
	fprintf(stderr, "wr = %i\n", trackPar->wR);
	fprintf(stderr, "wa = %i\n", trackPar->wA);
	fprintf(stderr, "navgR = %i\n", trackPar->navgR);
	fprintf(stderr, "navgA = %i\n", trackPar->navgA);
	fprintf(stderr, "edgepad = %i\n", trackPar->edgePad);
}

/*
  Parse up to 3 values from a RDF command list given a keyword using
  a format string. Number of values to parse is specified by n.
  Execution halts with appropriate error message if command not found
  or the data are invalid.
*/
static void parseValues(RDF *commands, char *keyword, char *format,
						void *d1, void *d2, void *d3, int32_t n)
{
	char *value;
	value = rdfValue(commands, keyword);
	if (value == NULL)
		error("%s not specified", keyword);
	if (sscanf(rdfValue(commands, keyword), format, d1, d2, d3) != n)
		error("Invalid %s", keyword);
}
