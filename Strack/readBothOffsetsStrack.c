#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "strack.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"
#include <sys/types.h>
#include "time.h"
#include <libgen.h>
#include "cRecipes/nrutil.h"

static void mallocOffsets(Offsets *offsets) {
	/*  Malloc buffers	*/
	int32_t i;
	float *fBuf, *fBufR;
	offsets->da = (float **)malloc(sizeof(float *) * offsets->na);
	offsets->dr = (float **)malloc(sizeof(float *) * offsets->na);
	fBuf = (float *)malloc(sizeof(float) * offsets->na * offsets->nr);
	fBufR = (float *)malloc(sizeof(float) * offsets->na * offsets->nr);
	for (i = 0; i < offsets->na; i++)
	{
		offsets->da[i] = &(fBuf[i * offsets->nr]);
		offsets->dr[i] = &(fBufR[i * offsets->nr]);
	}
}


static void readVrtMeta(GDALDatasetH hDS, Offsets *offsets)
{
	dictNode *metaData = NULL;
	float tmp;
	// Get meta data
	readDataSetMetaData(hDS, &metaData);
	// Write to offsets
	offsets->nr = GDALGetRasterXSize(hDS);
	offsets->na = GDALGetRasterYSize(hDS);
	// fprintf(stderr, "%s %s\n", get_value(metaData, "r0"), get_value(metaData, "a0"));
	offsets->rO = atoi(get_value(metaData, "r0"));
	offsets->aO = atoi(get_value(metaData, "a0"));
	// fprintf(stderr, "R0,A0\n");
	offsets->deltaR = atof(get_value(metaData, "deltaR"));
	offsets->deltaA = atof(get_value(metaData, "deltaA"));
}


int readBothOffsetsStrackVrt(Offsets *offsets, char *offsetsRoot) {
	char vrtBuffer[2048];
	GDALDatasetH hDS;
	int32_t status;
	const char *description;
	GDALRasterBandH hBand;
	
	if(checkForVrt(offsetsRoot, vrtBuffer) == NULL) return FALSE;
	hDS = GDALOpen(vrtBuffer, GDAL_OF_READONLY);
	readVrtMeta(hDS, offsets);
	mallocOffsets(offsets);
	
	// read range offsets
	hBand = GDALGetRasterBand(hDS, 1);
	description = GDALGetMetadataItem(hBand, "Description", NULL);
	
	fprintf(stderr, "Description %s\n", description);
	status = GDALRasterIO(hBand, GF_Read, 0, 0, offsets->nr, offsets->na, offsets->dr[0],
						  offsets->nr, offsets->na, GDT_Float32, 0, 0);
	if(status != 0 || strstr(description, "Range") == NULL) error("Problem reading range offsets from %s", vrtBuffer);
	
	// Read azimuth Offsets
	hBand = GDALGetRasterBand(hDS, 2);
	description = GDALGetMetadataItem(hBand, "Description", NULL);
	fprintf(stderr, "Description %s\n", description);
	status = GDALRasterIO(hBand, GF_Read, 0, 0, offsets->nr, offsets->na, offsets->da[0],
						  offsets->nr, offsets->na, GDT_Float32, 0, 0);
	if(status != 0 || strstr(description, "Azimuth") == NULL) error("Problem reading azimuth offsets from %s", vrtBuffer);
	return TRUE;

}


void readBothOffsetsStrack(Offsets *offsets)
{
	char *datFile, buf[1024], line[1024];
	FILE *fp;
	int32_t rO, aO, nr, na;
	float deltaA, deltaR;
	double c1;
	float *fBuf, *fBufR;
	char *file;
	int32_t lineCount, eod, i;
	/*  Read inputfile	*/
	file = offsets->file;
	datFile = &(buf[0]);
	buf[0] = '\0';
	datFile = strcat(datFile, file);
	datFile = strcat(datFile, ".dat");
	fp = openInputFile(datFile);
	lineCount = 0;
	lineCount = getDataString(fp, lineCount, line, &eod);
	if (sscanf(line, "%i%i%i%i%f%f", &rO, &aO, &nr, &na, &deltaR, &deltaA) != 6)
		error("%s  %i of %s", "readBothOffsetsStrack -- Missing image params at line:", lineCount, datFile);
	fclose(fp);
	/* Load offsets structure */
	offsets->nr = nr;
	offsets->na = na;
	offsets->rO = rO;
	offsets->aO = aO;
	offsets->deltaA = deltaA;
	offsets->deltaR = deltaR;
	
	mallocOffsets(offsets);
	/*  Read files	*/
	fprintf(stderr, "--- Reading offset file %s\n\n", file);
	fp = openInputFile(offsets->file);
	for (i = 0; i < na; i++)
		freadBS(offsets->da[i], sizeof(float), nr, fp, FLOAT32FLAG);
	fclose(fp);
	fprintf(stderr, "--- Reading offset file %s\n\n", offsets->rFile);
	fp = openInputFile(offsets->rFile);
	for (i = 0; i < na; i++)
		freadBS(offsets->dr[i], sizeof(float), nr, fp, FLOAT32FLAG);
	fclose(fp);
}