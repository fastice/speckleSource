#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "cullst.h"
#include "math.h"
#include <stdlib.h>
#include "mosaicSource/common/common.h"

#define STR_BUFFER_SIZE 1024
#define STR_BUFF(fmt, ...) ({                                    \
    char *__buf = (char *)calloc(STR_BUFFER_SIZE, sizeof(char)); \
    snprintf(__buf, STR_BUFFER_SIZE, fmt, ##__VA_ARGS__);        \
    __buf;                                                       \
})

static char **mallocByteMat(int32_t nA, int32_t nR);
static float **mallocFloatMat(int32_t nA, int32_t nR);
static void loadFile(void **data, char *filename, int32_t size, int32_t flag, int32_t nr, int32_t na);

static int32_t fileExists(char *vrtFile)
{
	if (access(vrtFile, F_OK) == 0)
	{
		return TRUE;
	}
	return FALSE;
}

void mallocBuffers(CullParams *cullPar) {
//	  space
	cullPar->offR = mallocFloatMat(cullPar->nA, cullPar->nR);
	cullPar->offA = mallocFloatMat(cullPar->nA, cullPar->nR);
	cullPar->offRS = mallocFloatMat(cullPar->nA, cullPar->nR);
	cullPar->offAS = mallocFloatMat(cullPar->nA, cullPar->nR);
	cullPar->corr = mallocFloatMat(cullPar->nA, cullPar->nR);
	cullPar->sigmaR = mallocFloatMat(cullPar->nA,cullPar->nR);
	cullPar->sigmaA = mallocFloatMat(cullPar->nA, cullPar->nR);
	cullPar->type = mallocByteMat(cullPar->nA, cullPar->nR);
}

void loadFromVRT(CullParams *cullPar) {
	GDALDatasetH hDS;
	dictNode *metaData = NULL, *metaDataMT=NULL;
	GDALRasterBandH hBand;
	int32_t nBands = 4, i;
	int32_t status;
	// Open data set
	hDS = GDALOpen(cullPar->inFileVRT, GDAL_OF_READONLY);
	// Get Size
	cullPar->nR = GDALGetRasterXSize(hDS);
	cullPar->nA = GDALGetRasterYSize(hDS);
	// Read and Process metadata
	readDataSetMetaData(hDS, &metaData);
	cullPar->metaData = metaData;
	cullPar->rStart = atoi(get_value(metaData, "r0"));
	cullPar->aStart = atoi(get_value(metaData, "a0"));
	cullPar->deltaR = atof(get_value(metaData, "deltaR"));
	cullPar->deltaA = atof(get_value(metaData, "deltaA"));
	fprintf(stderr, "\nParameters loaded\n");
	// Malloc buffers
	mallocBuffers(cullPar);
	fprintf(stderr, "Buffers Malloced\n");
	// Setup data for reading
	nBands = 3;
	int bands[] = {1, 2, 3};
	int dataTypes[] = {GDT_Float32, GDT_Float32, GDT_Float32};
	void *buffers[] = {cullPar->offR[0], cullPar->offA[0], cullPar->corr[0]};
	// Loop over and read bands
	for(i=0; i < nBands; i++) {
		hBand = GDALGetRasterBand(hDS, bands[i]);
		status = GDALRasterIO(hBand, GF_Read, 0, 0, cullPar->nR, cullPar->nA, buffers[i],
						 	  cullPar->nR, cullPar->nA, dataTypes[i], 0, 0);
		if(status != 0) error("loadFromVRT error reading band %i", i);
	}
	GDALClose(hDS);
	// Now match type
	GDALDatasetH hDSmt = GDALOpen(cullPar->inFileVRTMT, GDAL_OF_READONLY);
	readDataSetMetaData(hDSmt, &metaDataMT);
	cullPar->metaDataMT = metaDataMT;
	// Get MT Data
	hBand = GDALGetRasterBand(hDSmt, 1);
	status = GDALRasterIO(hBand, GF_Read, 0, 0, cullPar->nR, cullPar->nA, cullPar->type[0],
						 	  cullPar->nR, cullPar->nA, GDT_Byte, 0, 0);
	GDALClose(hDSmt);
}

static void populateMetaData(CullParams *cullPar)
{
	dictNode *metaData= NULL, *metaDataMT, *tmp;
	char filename[2048], geo2[2048];
	insert_node(&metaData, "r0", STR_BUFF("%i", cullPar->rStart));
    insert_node(&metaData, "a0", STR_BUFF("%i", cullPar->aStart ));
    insert_node(&metaData, "deltaA", STR_BUFF("%i", cullPar->deltaA));
    insert_node(&metaData, "deltaR", STR_BUFF("%i", cullPar->deltaR));
    insert_node(&metaData, "sigmaStreaks", STR_BUFF("%f", 0.));
    insert_node(&metaData, "sigmaRange", STR_BUFF("%f", 0.));
	if(cullPar->geo2 != NULL)
    	insert_node(&metaData, "geo1", STR_BUFF("%s", cullPar->geo1));
	if(cullPar->geo2 != NULL)
    	insert_node(&metaData, "geo2", STR_BUFF("%s", cullPar->geo2));
	for(tmp = metaData; tmp != NULL; tmp = tmp->next) 
		insert_node(&metaDataMT, tmp->key, STR_BUFF("%s", tmp->value));
	cullPar->metaDataMT = metaDataMT;
	// Assume MSB if loading using dat
	insert_node(&metaData, "ByteOrder", "MSB");
	cullPar->metaData = metaData;
}

static char **splitString(char *string) {
	char **words = malloc(sizeof(char *)*2);
	char *token = strtok(string, " ");
	words[0] = NULL; words[1] = NULL;
	int i = 0;
	while(token != NULL) {
		words[i] = token;
		token = strtok(NULL, " ");
		i++;
	}
	//Remove whitespace
	for( int j=0; j < i; j++) {
		words[j] = strtok(words[j], " \t\n");
	}
	return words;

}
void loadCullData(CullParams *cullPar)
{
	FILE *fp;
	char line[1024];
	char **split;
	int32_t lineCount, eod, i, j;
	int32_t r0, a0, nr, na;
	int32_t deltaA, deltaR;
	/*
	  Read .dat file
	*/
	fprintf(stderr, "Entering loadcull ...\n");
	if(fileExists(cullPar->inFileVRT)) {
		fprintf(stderr, "Found VRT, reading that ...\n");
		loadFromVRT(cullPar);
		return;
	}
	fprintf(stderr, "No VRT, using data files ...\n");
	fp = openInputFile(cullPar->inFileD);
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%i%i%i%i%i%i", &r0, &a0, &nr, &na, &deltaR, &deltaA) != 6)
		error("%s  %i of %s\n%s",
			  "readOffsets -- Missing image parameters at line:",
			  lineCount, cullPar->inFileD, line);
	lineCount = getDataString(fp, lineCount, line, &eod);
	if(strstr(line, "geodat") != NULL) {
		 split = splitString(line);
		 cullPar->geo1 = split[0];
		 cullPar->geo2 = split[1];
		 if(cullPar->geo1 != NULL && cullPar->geo2 !=NULL)
		 	fprintf(stderr, "%s %s\n", cullPar->geo1, cullPar->geo2);
	}
	fclose(fp);
	cullPar->nR = nr;
	cullPar->nA = na;
	cullPar->rStart = r0;
	cullPar->aStart = a0;
	cullPar->deltaA = deltaA;
	cullPar->deltaR = deltaR;
	//fprintf(stderr, "Populate Meta data...\n");
	populateMetaData(cullPar);
	//fprintf(stderr, "Print dict...\n");
	printDictionary(cullPar->metaData);
	//fprintf(stderr, "Malloc buffers...\n");
	mallocBuffers(cullPar);
	//fprintf(stderr, "Read dadta...\n");
	/*
	  Load data
	*/
	fprintf(stderr, "loading %s %i %i",cullPar->inFileR, nr, na );
	loadFile((void **)cullPar->offR, cullPar->inFileR, sizeof(float), nr, na, FLOAT32FLAG);
	fprintf(stderr, "loading %s %i %i",cullPar->inFileA, nr, na );
	loadFile((void **)cullPar->offA, cullPar->inFileA, sizeof(float), nr, na, FLOAT32FLAG);
	fprintf(stderr, "loading %s %i %i",cullPar->inFileC, nr, na );
	loadFile((void **)cullPar->corr, cullPar->inFileC, sizeof(float), nr, na, FLOAT32FLAG);
	fprintf(stderr, "loading %s %i %i",cullPar->inFileT, nr, na );
	loadFile((void **)cullPar->type, cullPar->inFileT, sizeof(char), nr, na, BYTEFLAG);
	for (i = 0; i < nr; i++)
		for (j = 0; j < na; j++)
			if (cullPar->offA[j][i] < (-LARGEINT + 1))
				cullPar->type[j][i] = 0;
	fprintf(stderr, "Leaving Loadcull....");
}

static void loadFile(void **data, char *filename, int32_t size, int32_t nr, int32_t na, int32_t flag)
{
	FILE *fp;
	fp = openInputFile(filename);
	freadBS(data[0], size, nr * na, fp, flag);
	fclose(fp);
}

void loadSimData(CullParams *cullPar)
{
	FILE *fp;
	char line[1024];
	int32_t lineCount, eod, i, j;
	int32_t r0, a0, nr, na;
	int32_t deltaA, deltaR;

	fp = openInputFile(cullPar->inFileD);
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%i%i%i%i%i%i", &r0, &a0, &nr, &na, &deltaR, &deltaA) != 6)
		error("%s  %i of %s\n%s", "readOffsets -- Missing image parameters at line:", lineCount, cullPar->inFileD, line);
	if (na != cullPar->nA || nr != cullPar->nR || r0 != cullPar->rStart || a0 != cullPar->aStart || deltaR != cullPar->deltaR || deltaA != cullPar->deltaA)
		error("UseSim only compatable if sim is same size as offset file");
	fclose(fp);
	fprintf(stderr, "Loading sim offsets\n");
	cullPar->offSimR = mallocFloatMat(na, nr);
	cullPar->offSimA = mallocFloatMat(na, nr);
	/* Read offsets */
	loadFile((void **)cullPar->offSimR, cullPar->inFileSimR, sizeof(float), nr, na, FLOAT32FLAG);
	loadFile((void **)cullPar->offSimA, cullPar->inFileSimA, sizeof(float), nr, na, FLOAT32FLAG);
}

int32_t loadCullMask(CullParams *cullPar)
{
	FILE *fp, *fpD;
	char line[1024];
	int32_t lineCount, eod, i, j;
	int32_t r0, a0, nr, na;
	int32_t deltaA, deltaR;

	fprintf(stderr, "Found offsets.mask\n");
	fprintf(stderr, "%d %d\n", cullPar->nA, cullPar->nR);
	/* Read offsets.dat and check the file is the same - mostly to avoid mask register offsets - added 08/23/18 */
	fpD = openInputFile("offsets.dat");
	lineCount = getDataString(fpD, lineCount, line, &eod);
	if (sscanf(line, "%i%i%i%i%i%i", &r0, &a0, &nr, &na, &deltaR, &deltaA) != 6)
		error("%s  %i of %s\n%s", "readOffsets -- Missing image parameters at line:", lineCount, "offsets.dat", line);
	fclose(fpD);

	if (cullPar->nA != na || cullPar->nR != nr)
	{
		fprintf(stderr, "Ignoring mask because it is not the same size: regoffsets ?\n");
		return (FALSE);
	}
	cullPar->mask = mallocByteMat(cullPar->nA, cullPar->nR);
	loadFile((void **)cullPar->mask, "offsets.mask", sizeof(char), nr, na, BYTEFLAG);
	return TRUE;
}

static char **mallocByteMat(int32_t nA, int32_t nR)
{
	char *tmp, **tmp1;

	int32_t i;
	tmp = calloc(nR * nA, sizeof(char));
	tmp1 = (char **)calloc(nA, sizeof(char *));
	for (i = 0; i < nA; i++)
		tmp1[i] = &(tmp[i * nR]);
	return tmp1;
}

static float **mallocFloatMat(int32_t nA, int32_t nR)
{
	float *tmp, **tmp1;
	int32_t i;
	size_t size;
	size = (size_t)nR * (size_t)nA;
	tmp = calloc(size, sizeof(float));
	tmp1 = (float **)calloc(nA, sizeof(float *));
	for (i = 0; i < nA; i++)
		tmp1[i] = &(tmp[i * nR]);
	return tmp1;
}
