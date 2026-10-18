#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "clib/standard.h"
#include "strack.h"
#include "math.h"
#include "sys/stat.h"
#include "gdal.h"
#include "cpl_conv.h"
#include "gdalIO/gdalIO/grimpgdal.h"

/*  Read the interferogram through GDAL when a VRT backs it.  Returns TRUE if it handled the
    file.  This is the NISAR path: a RIFG converted to complex float32 and cropped to the SLC
    subset, with everything getInt otherwise takes from a legacy geodat carried as VRT
    metadata (nisarhdf writes these keys):

        MLRangeSize MLAzimuthSize NumberRangeLooks NumberAzimuthLooks Wavelength centerLat

    plus, written by the cropping step, the SLC pixel of interferogram sample (0,0):

        SLCRangeStart SLCAzimuthStart        (absent -> 0, i.e. an uncropped interferogram)

    The legacy complex-int16 + geodat path below is untouched. */
static int32_t getIntVrt(TrackParams *trackPar, char *vrtFile)
{
	GDALDatasetH hDS;
	GDALRasterBandH hBand;
	char **meta;
	const char *v;
	fftw_complex *buf;
	float *line;
	int32_t i, j, nr, na;

	hDS = GDALOpen(vrtFile, GDAL_OF_READONLY);
	if (hDS == NULL)
	{
		return (FALSE);
	}
	hBand = GDALGetRasterBand(hDS, 1);
	if (GDALGetRasterDataType(hBand) != GDT_CFloat32)
	{
		fprintf(stderr, "getIntVrt: %s is not CFloat32, falling back to the legacy reader\n", vrtFile);
		GDALClose(hDS);
		return (FALSE);
	}
	nr = GDALGetRasterBandXSize(hBand);
	na = GDALGetRasterBandYSize(hBand);
	meta = GDALGetMetadata(hDS, NULL);
	trackPar->intDat.nr = nr;
	trackPar->intDat.na = na;
	trackPar->intDat.nrl = 1;
	trackPar->intDat.nal = 1;
	trackPar->intDat.r0 = 0;
	trackPar->intDat.a0 = 0;
	if ((v = CSLFetchNameValue(meta, "NumberRangeLooks")) != NULL)
	{
		trackPar->intDat.nrl = atoi(v);
	}
	if ((v = CSLFetchNameValue(meta, "NumberAzimuthLooks")) != NULL)
	{
		trackPar->intDat.nal = atoi(v);
	}
	if ((v = CSLFetchNameValue(meta, "SLCRangeStart")) != NULL)
	{
		trackPar->intDat.r0 = atoi(v);
	}
	if ((v = CSLFetchNameValue(meta, "SLCAzimuthStart")) != NULL)
	{
		trackPar->intDat.a0 = atoi(v);
	}
	if ((v = CSLFetchNameValue(meta, "Wavelength")) != NULL)
	{
		trackPar->lambda = atof(v);
	}
	if ((v = CSLFetchNameValue(meta, "centerLat")) != NULL)
	{
		trackPar->latc = atof(v);
	}
	if (trackPar->intDat.nrl < 1 || trackPar->intDat.nal < 1)
	{
		error("getIntVrt: bad looks in %s", vrtFile);
	}
	buf = (fftw_complex *)strackMallocComplex((size_t)nr * (size_t)na);
	if (buf == NULL)
	{
		fprintf(stderr, "getIntVrt: cannot malloc interferogram buffer\n");
		GDALClose(hDS);
		trackPar->intDat.intf = NULL;
		return (TRUE);
	}
	trackPar->intDat.intf = (fftw_complex **)malloc((size_t)na * sizeof(fftw_complex *));
	line = (float *)malloc((size_t)nr * 2 * sizeof(float));
	for (i = 0; i < na; i++)
	{
		trackPar->intDat.intf[i] = &(buf[(size_t)i * (size_t)nr]);
		if (GDALRasterIO(hBand, GF_Read, 0, i, nr, 1, line, nr, 1, GDT_CFloat32, 0, 0) != CE_None)
		{
			error("getIntVrt: read failed on %s", vrtFile);
		}
		for (j = 0; j < nr; j++)
		{
			trackPar->intDat.intf[i][j].re = line[2 * j];
			trackPar->intDat.intf[i][j].im = line[2 * j + 1];
		}
	}
	free(line);
	GDALClose(hDS);
	fprintf(stderr, "getIntVrt: %s  nr %i na %i looks %ix%i  SLC origin r0 %i a0 %i  lambda %f\n",
			vrtFile, nr, na, trackPar->intDat.nrl, trackPar->intDat.nal,
			trackPar->intDat.r0, trackPar->intDat.a0, trackPar->lambda);
	return (TRUE);
}
/*
  Malloc space, read params for, and input complex interferogram
*/

void getInt(TrackParams *trackPar)
{
	inputImageStructure inputImage;
	fftw_complex *buf;
	size_t bufSize;
	strackComplex *inBuf;
	FILE *fp;
	struct stat st;
	int32_t i, j, nLines;

	char vrtBuf[4096], *vrtFile;

	fprintf(stderr, "get int \n");
	/*  VRT-backed (NISAR RIFG) interferogram: everything comes from the file and its metadata,
	    so no geodat is needed.  Falls through to the legacy reader if there is no VRT. */
	vrtFile = checkForVrt(trackPar->intFile, vrtBuf);
	if (vrtFile != NULL)
	{
		if (getIntVrt(trackPar, vrtFile) == TRUE)
		{
			return;
		}
	}
	if (trackPar->intGeodat == NULL)
		error("getInt: Missing geodat filename");

	inputImage.stateFlag = TRUE;
	parseInputFile(trackPar->intGeodat, &inputImage);
	trackPar->latc = inputImage.latControlPoints[0];
	trackPar->intDat.nr = inputImage.rangeSize;
	trackPar->intDat.na = inputImage.azimuthSize;
	trackPar->intDat.nrl = inputImage.nRangeLooks;
	trackPar->intDat.nal = inputImage.nAzimuthLooks;
	trackPar->intDat.r0 = 0;
	trackPar->intDat.a0 = 0;
	trackPar->lambda = inputImage.par.lambda;
	inBuf = (strackComplex *)malloc(trackPar->intDat.nr * sizeof(strackComplex));
	trackPar->intDat.intf = (fftw_complex **)malloc(trackPar->intDat.na * sizeof(fftw_complex *));
	bufSize = (int32_t)trackPar->intDat.na * (int32_t)trackPar->intDat.nr * (int32_t)sizeof(fftw_complex);
	buf = (fftw_complex *)strackMallocComplex(bufSize / sizeof(fftw_complex));
	if (buf == NULL)
	{
		fprintf(stderr, "Cannot malloc interferogram buffer, proceeding with no interferogram\n");
		trackPar->intDat.intf = NULL;
		return;
	}
	fprintf(stderr, "Reading Interferogram nr,na, bufsize %d %d %zu \n", trackPar->intDat.nr, trackPar->intDat.na, bufSize);

	for (i = 0; i < trackPar->intDat.na; i++)
	{
		trackPar->intDat.intf[i] = &(buf[i * trackPar->intDat.nr]);
		/* Fill buffer with something close to zero */
		for (j = 0; j < trackPar->intDat.nr; j++)
		{
			trackPar->intDat.intf[i][j].re = 1;
			trackPar->intDat.intf[i][j].im = 1;
		}
	}

	/* Allow interferogram to be smaller than image - only read in part with data */
	stat(trackPar->intFile, &st);
	/* Actual number of valid lines */
	nLines = st.st_size / (trackPar->intDat.nr * 4);
	fprintf(stderr, "geodat lines %i int %i\n", trackPar->intDat.na, nLines);
	fp = fopen(trackPar->intFile, "r");
	/*
	  Read interferogram - if short the last part will have fill values defined above.
	 */
	for (i = 0; i < nLines; i++)
	{
		freadBS(inBuf, sizeof(strackComplex), trackPar->intDat.nr, fp, INT16FLAG);
		for (j = 0; j < trackPar->intDat.nr; j++)
		{
			trackPar->intDat.intf[i][j].re = (float)inBuf[j].r;
			trackPar->intDat.intf[i][j].im = (float)inBuf[j].i;
		}
	}
	fclose(fp);
	free(inBuf);
}
