#include <stdio.h>
#include "string.h"
#include "clib/standard.h"
#include "math.h"
#include <stdlib.h>
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "cullls.h"
#include "gdalIO/gdalIO/grimpgdal.h"
static void readArgs(int argc, char *argv[], CullLSParams *cullPar);
static void usage();
/*
   Global variables definitions
*/
int32_t RangeSize = 0;		 /* Range size of complex image */
int32_t AzimuthSize = 0;	 /* Azimuth size of complex image */
int32_t BufferSize = 0;		 /* Size of nonoverlap region of the buffer */
int32_t BufferLines = 0;	 /* # of lines of nonoverlap in buffer */
double RangePixelSize = 0;	 /* Range PixelSize */
double AzimuthPixelSize = 0; /* Azimuth PixelSize */
double SLat = -91.;
static void writeLSCulledOffsets(CullLSParams *cullPar);
static void writeLSCulledTiffs(CullLSParams *cullPar, int32_t nx, int32_t ny);

/*
 * gdalIO.o's checkForVrt() references appendSuffix, which lives in
 * mosaicSource/common/initRoutines.c. cullls does not link $(COMMON), so it
 * needs its own definition to satisfy the linker, exactly as
 * utilsSource/intFloat/vrtIO.c and gdalIO/test/testio.c already do. cullls
 * never calls checkForVrt itself.
 */
char *appendSuffix(char *file, char *suffix, char *buf)
{
	strcpy(buf, file);
	strcat(buf, suffix);
	return buf;
}

int main(int argc, char *argv[])
{
	CullLSParams cullPar;

	/* Every other GDAL-using program in the tree registers the drivers at the
	   top of main (strack.c:59, cullst.c:39, ...). cullls did not, because it
	   had no GDAL output until -tiff. */
	GDALAllRegister();
	readArgs(argc, argv, &cullPar);
	fprintf(stderr, "cullIslandThresh %i\n", cullPar.islandThresh);
	/*
	   Load data
	*/
	loadLSCullData(&cullPar);
	/*
	  cull data
	*/
	fprintf(stderr, "Cull data\n");
	cullPar.nMatch = -1; /* Use negative value so this only gets set on the first pass through */
	cullLSData(&cullPar);
	cullLSData(&cullPar);
	cullLSData(&cullPar);
	/*
	  Compute stats for data
	*/
	fprintf(stderr, "Cull Stats\n");
	cullLSStats(&cullPar);
	/*
	  Compute stats for data
	*/
	fprintf(stderr, "Cull Smooth\n");
	cullLSSmooth(&cullPar);
	/*
	  cull islands
	*/
	if (cullPar.islandThresh > 0)
		cullLSIslands(&cullPar);
	/*
	  Output result
	*/
	fprintf(stderr, "Cull write\n");
	writeLSCulledOffsets(&cullPar);
}


/******************************************************************************************
   Write the culled products as GeoTIFFs plus a wrapping VRT.

   Unlike the SAR cull products, Landsat offsets are genuinely map projected: the
   .dat carries x0/y0 in metres, the base pixel size dx/dy and the match step, and
   an EPSG code. So these get a real north-up geotransform and a projection via
   saveAsGeotiff, not the pixel-coordinate geotransform writeFlatTiff stamps on
   radar-geometry rasters.

   Names simply append .tif to the raw name (match.X.Y.cull.dx.tif), which is the
   convention every Python writer in the chain uses and what intfloat -inputVRT
   and the mosaic reader expect. That is deliberately not Cullst's deriveTif(),
   which strips the extension and substitutes a band role.
*******************************************************************************************/
static void writeLSCulledTiffs(CullLSParams *cullPar, int32_t nx, int32_t ny)
{
	double geoTransform[6];
	char epsg[32], tifName[2048], vrtName[2048];
	const char *bands[4];
	const char *bandNames[4] = {"dx", "dy", "sx", "sy"};
	float noData[4] = {(float)NODATA, (float)NODATA, (float)NODATA, (float)NODATA};
	char dxTif[2048], dyTif[2048], sxTif[2048], syTif[2048];
	/* The match posting is the base pixel size times the tracking step */
	double deltaX = cullPar->matches.dx * (double)cullPar->matches.stepX;
	double deltaY = cullPar->matches.dy * (double)cullPar->matches.stepY;
	/*
	  x0/y0 are the lower-left pixel centre in metres, which is what
	  computeGeoTransform expects; saveAsGeotiff flips to north-up on write.
	*/
	computeGeoTransform(geoTransform, cullPar->matches.x0, cullPar->matches.y0,
						nx, ny, deltaX, deltaY);
	sprintf(epsg, "%i", cullPar->fitDat.proj);
	/*
	  No dataset metadata: the .cull.dat sidecar is still written alongside and
	  carries the acquisition dates, sigmas and rates that downstream code reads.
	*/
	/*
	  One GeoTIFF per band, then a VRT over the four float bands. The match-type
	  byte band is written separately so the float VRT stays homogeneous.
	*/
	sprintf(dxTif, "%s.cull.dx.tif", cullPar->fitDat.matchFile);
	sprintf(dyTif, "%s.cull.dy.tif", cullPar->fitDat.matchFile);
	sprintf(sxTif, "%s.cull.sx.tif", cullPar->fitDat.matchFile);
	sprintf(syTif, "%s.cull.sy.tif", cullPar->fitDat.matchFile);
	saveAsGeotiff(dxTif, cullPar->XS[0], nx, ny, geoTransform, epsg, NULL,
				  "GTiff", GDT_Float32, (float)NODATA);
	saveAsGeotiff(dyTif, cullPar->YS[0], nx, ny, geoTransform, epsg, NULL,
				  "GTiff", GDT_Float32, (float)NODATA);
	saveAsGeotiff(sxTif, cullPar->matches.sigmaX[0], nx, ny, geoTransform, epsg,
				  NULL, "GTiff", GDT_Float32, (float)NODATA);
	saveAsGeotiff(syTif, cullPar->matches.sigmaY[0], nx, ny, geoTransform, epsg,
				  NULL, "GTiff", GDT_Float32, (float)NODATA);
	sprintf(tifName, "%s.cull.mtype.tif", cullPar->fitDat.matchFile);
	saveAsGeotiff(tifName, cullPar->matches.type[0], nx, ny, geoTransform, epsg,
				  NULL, "GTiff", GDT_Byte, 0.);
	bands[0] = dxTif; bands[1] = dyTif; bands[2] = sxTif; bands[3] = syTif;
	sprintf(vrtName, "%s.cull.vrt", cullPar->fitDat.matchFile);
	makeTiffVRTNamed(vrtName, bands, bandNames, 4, noData, NULL);
	fprintf(stderr, "%s\n", vrtName);
}

/******************************************************************************************
   Write culled offsets
*******************************************************************************************/
static void writeLSCulledOffsets(CullLSParams *cullPar)
{
	uint32_t i, j;
	FILE *fp;
	char *file1;
	size_t sl;
	int32_t nx, ny;
	extern int32_t nMatch, nAttempt, nTotal;
	int32_t nGoodPts;
	double sigmaXAvg, sigmaYAvg;
	sl = 1500;
	file1 = (char *)malloc(sl);

	nx = cullPar->matches.nx;
	ny = cullPar->matches.ny;
	nGoodPts = 0;
	/*
	  Final count
	*/
	sigmaXAvg = 0;
	sigmaYAvg = 0;
	for (i = 0; i < ny; i++)
	{
		for (j = 0; j < nx; j++)
		{
			if (cullPar->XS[i][j] > (NODATA + 1) && cullPar->YS[i][j] > (NODATA + 1))
			{
				sigmaXAvg += cullPar->matches.sigmaX[i][j] * cullPar->matches.sigmaX[i][j];
				sigmaYAvg += cullPar->matches.sigmaY[i][j] * cullPar->matches.sigmaY[i][j];
				nGoodPts++;
			}
		}
	}
	sigmaXAvg = sqrt(sigmaXAvg / (double)max(nGoodPts, 1));
	sigmaYAvg = sqrt(sigmaYAvg / (double)max(nGoodPts, 1));
	/*	fprintf(stderr,"File root %s %i %i\n",matchP->outputFile, nx,ny);
		for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,matchP->outputFile); file1=strcat(file1,".rho"); fprintf(stderr,"%s\n",file1);
		fp=fopen(file1,"w");  fwriteBS(matches->Rho[0],sizeof(float),(size_t)(nx*ny),fp,FLOAT32FLAG); fclose(fp);*/
	if (cullPar->tiffFlag == TRUE)
	{
		writeLSCulledTiffs(cullPar, nx, ny);
		goto writeDat;
	}
	for (i = 0; i < sl; i++)
		file1[i] = '\0';
	file1 = strcpy(file1, cullPar->fitDat.matchFile);
	file1 = strcat(file1, ".cull.dx");
	fprintf(stderr, "%s\n", file1);
	fp = fopen(file1, "w");
	fwriteBS(cullPar->XS[0], sizeof(float), (size_t)(nx * ny), fp, FLOAT32FLAG);
	fclose(fp);

	for (i = 0; i < sl; i++)
		file1[i] = '\0';
	file1 = strcpy(file1, cullPar->fitDat.matchFile);
	file1 = strcat(file1, ".cull.dy");
	fprintf(stderr, "%s\n", file1);
	fp = fopen(file1, "w");
	fwriteBS(cullPar->YS[0], sizeof(float), (size_t)(nx * ny), fp, FLOAT32FLAG);
	fclose(fp);

	for (i = 0; i < sl; i++)
		file1[i] = '\0';
	file1 = strcpy(file1, cullPar->fitDat.matchFile);
	file1 = strcat(file1, ".cull.sx");
	fprintf(stderr, "%s\n", file1);
	fp = fopen(file1, "w");
	fwriteBS(cullPar->matches.sigmaX[0], sizeof(float), (size_t)(nx * ny), fp, FLOAT32FLAG);
	fclose(fp);

	for (i = 0; i < sl; i++)
		file1[i] = '\0';
	file1 = strcpy(file1, cullPar->fitDat.matchFile);
	file1 = strcat(file1, ".cull.sy");
	fprintf(stderr, "%s\n", file1);
	fp = fopen(file1, "w");
	fwriteBS(cullPar->matches.sigmaY[0], sizeof(float), (size_t)(nx * ny), fp, FLOAT32FLAG);
	fclose(fp);

	for (i = 0; i < sl; i++)
		file1[i] = '\0';
	file1 = strcpy(file1, cullPar->fitDat.matchFile);
	file1 = strcat(file1, ".cull.mtype");
	fprintf(stderr, "%s\n", file1);
	fp = fopen(file1, "w");
	fwrite(cullPar->matches.type[0], sizeof(uint8_t), (size_t)(nx * ny), fp);
	fclose(fp);

writeDat:
	/* The .dat sidecar is written in both modes: it carries the acquisition
	   dates, slowFlag, sigmas and match rates that a geotransform cannot hold
	   and that mosaic3d, lsdat and velscreen all read. */
	for (i = 0; i < sl; i++)
		file1[i] = '\0';
	file1 = strcpy(file1, cullPar->fitDat.matchFile);
	file1 = strcat(file1, ".cull.dat");
	fprintf(stderr, "%s\n", file1);

	fp = fopen(file1, "w");
	fprintf(fp, "fileEarly = %s\n", cullPar->matches.fileEarly);
	fprintf(fp, "fileLate = %s\n", cullPar->matches.fileLate);
	fprintf(fp, "x0 = %11.2lf\n", cullPar->matches.x0);
	fprintf(fp, "y0 = %11.2lf\n", cullPar->matches.y0);
	fprintf(fp, "dx = %10.5lf\n", cullPar->matches.dx);

	fprintf(fp, "dy = %10.5lf\n", cullPar->matches.dy);
	fprintf(fp, "stepY = %u\n", cullPar->matches.stepX);
	fprintf(fp, "stepX = %u\n", cullPar->matches.stepY);
	fprintf(fp, "nx = %9i\n", nx);
	fprintf(fp, "ny = %9i\n", ny);
	fprintf(fp, "slowFlag = %9i\n", cullPar->fitDat.slowFlag);
	fprintf(fp, "EPSG = %9i\n", cullPar->fitDat.proj);
	fprintf(fp, "earlyImageJD = %10lf\n", cullPar->matches.jdEarly);
	fprintf(fp, "lateImageJD = %10lf\n", cullPar->matches.jdLate);
	fprintf(fp, "IntervalBetweenImages = %5i\n", (int)(cullPar->matches.jdLate - cullPar->matches.jdEarly));
	fprintf(fp, "Success_rate_for_attempted_matches(%%) =  %7.2f \n", ((float)cullPar->nMatch / (float)max(cullPar->nAttempt, 1)) * 100.);
	fprintf(fp, "Culled_rate_for_attempted_matches(%%) =  %7.2f \n", ((float)nGoodPts / (float)max(cullPar->nAttempt, 1)) * 100.);
	fprintf(fp, "Mean_sigmaX = %11.3lf\n", sigmaXAvg);
	fprintf(fp, "Mean_sigmaY = %11.3lf\n", sigmaYAvg);
	fprintf(fp, "& \n");
	fclose(fp);
}

/******************************************************************************************
	 read args
*******************************************************************************************/
static void readArgs(int argc, char *argv[], CullLSParams *cullPar)
{
	int32_t filenameArg;
	char *argString;
	char *inBase, *outBase;
	int32_t islandThresh;
	int32_t sx, sy;
	int32_t boxSize, nGood;
	float maxY, maxX;
	int32_t i, n, sLen;

	if (argc < 2 || argc > 19)
		usage(); /* Check number of args */
	/* prog + trailing inbase; the loop consumes its own flag values, so this
	   also admits flags that take no value (-tiff). It was argc - 3, which
	   silently dropped a valueless final flag. */
	n = argc - 2;
	sx = 3;
	sy = 3;
	boxSize = 9;
	nGood = 17;
	maxY = 1.0;
	maxX = .75;
	islandThresh = -1;
	cullPar->tiffFlag = FALSE;
	for (i = 1; i <= n; i++)
	{
		argString = strchr(argv[i], '-');
		if (strstr(argString, "sx") != NULL)
		{
			sscanf(argv[i + 1], "%i", &sx);
			i++;
		}
		else if (strstr(argString, "sy") != NULL)
		{
			sscanf(argv[i + 1], "%i", &sy);
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
		else if (strstr(argString, "maxY") != NULL)
		{
			sscanf(argv[i + 1], "%f", &maxY);
			i++;
		}
		else if (strstr(argString, "maxX") != NULL)
		{
			sscanf(argv[i + 1], "%f", &maxX);
			i++;
		}
		else if (strstr(argString, "islandThresh") != NULL)
		{
			sscanf(argv[i + 1], "%i", &islandThresh);
			i++;
		}
		else if (strstr(argString, "tiff") != NULL)
		{
			cullPar->tiffFlag = TRUE;
		}
		else
			usage();
	}
	cullPar->islandThresh = islandThresh;
	cullPar->sX = sx;
	cullPar->sY = sy;
	cullPar->bR = boxSize;
	cullPar->bA = boxSize;
	cullPar->nGood = nGood;
	cullPar->maxY = maxY;
	cullPar->maxX = maxX;
	cullPar->fitDat.matchFile = argv[argc - 1];
	return;
}

static void usage()
{
	error("cullls -islandThresh islandThresh -maxX maxX -maxY maxY -nGood nGood -boxSize boxSize -sx sx -sy sy -tiff inbase\n%s\n\n%s\n",
		  "where",
		  "sx,sy = smoothing window size to apply",
		  "islandThresh = cull isolated islands < islandThresh in diameter",
		  "maxX,maxY = max deviation from local median",
		  "tiff = write culled products as GeoTIFF + VRT rather than raw flat files");
}
