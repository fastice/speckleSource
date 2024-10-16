#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "strack.h"
#include "math.h"
//#include "gdalIO/gdalIO/grimpgdal.h"
/*
  Malloc space, read params for, and input complex interferogram
*/


static void mallocMaskBuff(TrackParams *trackPar)
{
	int32_t i;
	char *buf;
	trackPar->maskDat.mask = (char **)malloc(trackPar->maskDat.na * sizeof(char *));
	buf = (char *)malloc(trackPar->maskDat.na * trackPar->maskDat.nr * sizeof(char));
	for (i = 0; i < trackPar->maskDat.na; i++)
		trackPar->maskDat.mask[i] = &(buf[i * trackPar->maskDat.nr]);
}

static void getMaskVRT(TrackParams *trackPar) {
	GDALDatasetH hDS;
	dictNode *metaData = NULL;
	GDALRasterBandH hBand;
	int32_t status;
	float tmp;
	fprintf(stderr, "Reading mask specified by vrt %s", trackPar->maskvrt);
	// Get meta data
	if (access(trackPar->maskvrt, F_OK) != 0) error("getMaskVRT: Cannot open vrt %s", trackPar->maskvrt);
	hDS = GDALOpen(trackPar->maskvrt, GDAL_OF_READONLY);
	hBand = GDALGetRasterBand(hDS, 1);
	readDataSetMetaData(hDS, &metaData);	
	trackPar->maskDat.nr = GDALGetRasterBandXSize(hBand);
	trackPar->maskDat.na  = GDALGetRasterBandYSize(hBand);
	trackPar->maskDat.r0 = atoi(get_value(metaData, "r0"));
	trackPar->maskDat.a0 = atoi(get_value(metaData, "a0"));
	// nrl is the skip, for offsets
	trackPar->maskDat.nrl = atoi(get_value(metaData, "deltaR"));
	trackPar->maskDat.nal = atoi(get_value(metaData, "deltaA"));
	// Malloc image space
	fprintf(stderr, "Mallocing space for mask...\n");
	mallocMaskBuff(trackPar);
	// Read the data
	fprintf(stderr, "Reading mask...\n");
	status = GDALRasterIO(hBand, GF_Read, 0, 0, trackPar->maskDat.nr , trackPar->maskDat.na, trackPar->maskDat.mask[0],
						  trackPar->maskDat.nr , trackPar->maskDat.na, GDT_Byte, 0, 0);
	fprintf(stderr, "Mask read\n");
}

void getMask(TrackParams *trackPar)
{
	inputImageStructure inputImage;
	char *buf, *file, *datFile, buf1[2048];
	FILE *fp;
	int32_t i, j;
	int32_t r0, a0, nr, na;
	float deltaR, deltaA;
	int32_t lineCount, eod;
	char line[1024];
	// If there is a vrt, use it.
	if(trackPar->maskvrt != NULL)
	{
		getMaskVRT(trackPar);
		return;
	}
	// Otherwise use the dat file
	if (trackPar->maskGeodat == NULL)
		error("getMask: Missing geodat filename");
	inputImage.stateFlag = TRUE;
	if (trackPar->maskType == GEODATMASK)
	{
		parseInputFile(trackPar->maskGeodat, &inputImage);
		trackPar->maskDat.r0 = 0;
		trackPar->maskDat.a0 = 0;
		trackPar->maskDat.nr = inputImage.rangeSize;
		trackPar->maskDat.na = inputImage.azimuthSize;
		trackPar->maskDat.nrl = inputImage.nRangeLooks;
		trackPar->maskDat.nal = inputImage.nAzimuthLooks;
	}
	else if (trackPar->maskType == OFFSETMASK)
	{
		file = trackPar->maskGeodat;
		fp = openInputFile(file);
		lineCount = 0;
		lineCount = getDataString(fp, lineCount, line, &eod);
		if (sscanf(line, "%i%i%i%i%f%f", &r0, &a0, &nr, &na, &deltaR, &deltaA) != 6)
			error("%s  %i of %s", "readOffsets -- Missing image parameters at line:", lineCount, datFile);
		fclose(fp);
		trackPar->maskDat.r0 = r0;
		trackPar->maskDat.a0 = a0;
		trackPar->maskDat.nr = nr;
		trackPar->maskDat.na = na;
		trackPar->maskDat.nrl = (int)deltaR;
		trackPar->maskDat.nal = (int)deltaA;
		
	}
	else
		error("Invalid mask type");
	fprintf(stderr, "\nmask params r0,a0,nr,na,deltaR,deltaA %i %i %i %i %i %i\n",
			trackPar->maskDat.r0, trackPar->maskDat.a0, trackPar->maskDat.nr, trackPar->maskDat.na,
			trackPar->maskDat.nrl, trackPar->maskDat.nal);
	// Allocate space for mask
	mallocMaskBuff(trackPar);
	// Open and read file
	fp = fopen(trackPar->maskFile, "r");
	fprintf(stderr, "Reading file %s %i %i\n", trackPar->maskFile, trackPar->maskDat.nr, trackPar->maskDat.na);
	freadBS(trackPar->maskDat.mask[0], sizeof(char), trackPar->maskDat.nr * trackPar->maskDat.na, fp, BYTEFLAG);

	for (i = 0; i < trackPar->maskDat.na; i += 200)
	{
		for (j = 0; j < trackPar->maskDat.nr; j += 100)
			fprintf(stderr, "%1i", trackPar->maskDat.mask[i][j]);
		fprintf(stderr, "   %10i\n", i * trackPar->maskDat.nal + trackPar->maskDat.a0);
	}

	fclose(fp);
}
