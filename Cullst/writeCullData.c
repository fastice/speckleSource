#include "stdio.h"
#include "string.h"
#include <stdlib.h>
#include "clib/standard.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "cullst.h"
#include "math.h"

#define LSB 0
#define MSB 1

// From writeVrt.c in strack
//void writeSingleVRT(int32_t nR, int32_t nA, dictNode *metaData, char *vrtFile, char *bandFiles[], char *bandNames[],
//                    GDALDataType dataTypes[], char *byteSwapOption, int32_t nBands);

/* Derive a per-band GeoTIFF name from a raw band filename: strip the trailing
   ".xx" extension and append ".<role>.tif" (so makeTiffVRT's filename-derived
   band Description becomes <role>, matching the raw-VRT band names). */
static char *deriveTif(const char *base, const char *role)
{
    char *out = (char *)malloc(strlen(base) + strlen(role) + 8);
    char *dot;
    strcpy(out, base);
    dot = strrchr(out, '.');
    if (dot != NULL)
        *dot = '\0';
    strcat(out, ".");
    strcat(out, role);
    strcat(out, ".tif");
    return out;
}

/* Write a flat row-major buffer to a GeoTIFF with a pixel-coord geotransform.
   No vertical flip: row 0 of the buffer becomes the top row of the tif, matching
   the raw-binary / writeSingleVRT convention (row 0 = azimuth 0). */
static void writeCullTiff(const char *filename, void *flatData, int32_t width,
                          int32_t height, GDALDataType dataType,
                          float noDataValue, dictNode *metaData)
{
    double pixGT[6] = {-0.5, 1., 0., -0.5, 0., 1.};
    const char *options[] = {"COMPRESS=DEFLATE", NULL};
    GDALDriverH driver = GDALGetDriverByName("GTiff");
    GDALDatasetH ds = GDALCreate(driver, filename, width, height, 1, dataType, (char **)options);
    if (ds == NULL) {
        fprintf(stderr, "writeCullTiff: cannot create %s\n", filename);
        return;
    }
    GDALSetGeoTransform(ds, pixGT);
    GDALRasterBandH band = GDALGetRasterBand(ds, 1);
    GDALSetRasterNoDataValue(band, noDataValue);
    GDALRasterIO(band, GF_Write, 0, 0, width, height, flatData, width, height, dataType, 0, 0);
    if (metaData != NULL) {
        writeDataSetMetaData(ds, metaData);
    }
    GDALClose(ds);
}

/* Write the cull outputs as per-band GeoTIFFs (Range/Azimuth/RangeSigma/
   AzimuthSigma/Correlation + MatchType), wrapped by the same .vrt / .mt.vrt the
   raw path writes. The sigma tiffs keep the raw path's inBase location. */
static void writeCullTiffData(CullParams *cullPar)
{
    int32_t nr = cullPar->nR, na = cullPar->nA;
    char *rTif = deriveTif(cullPar->outFileR, "RangeOffsets");
    char *aTif = deriveTif(cullPar->outFileA, "AzimuthOffsets");
    char *srTif = deriveTif(cullPar->outFileSR, "RangeSigma");
    char *saTif = deriveTif(cullPar->outFileSA, "AzimuthSigma");
    char *cTif = deriveTif(cullPar->outFileC, "Correlation");
    char *tTif = deriveTif(cullPar->outFileT, "MatchType");
    /* No vertical flip: row 0 = azimuth 0, matching the raw convention. */
    writeCullTiff(rTif, cullPar->offRS[0], nr, na, GDT_Float32, -2.e9f, cullPar->metaData);
    writeCullTiff(aTif, cullPar->offAS[0], nr, na, GDT_Float32, -2.e9f, cullPar->metaData);
    writeCullTiff(srTif, cullPar->sigmaR[0], nr, na, GDT_Float32, -2.e9f, cullPar->metaData);
    writeCullTiff(saTif, cullPar->sigmaA[0], nr, na, GDT_Float32, -2.e9f, cullPar->metaData);
    writeCullTiff(cTif, cullPar->corr[0], nr, na, GDT_Float32, -2.e9f, cullPar->metaData);
    writeCullTiff(tTif, cullPar->type[0], nr, na, GDT_Byte, 0.f, cullPar->metaDataMT);
    /* Wrap: main VRT = R/A/SR/SA/C (5 bands); .mt.vrt = MatchType. */
    const char *mainBands[5] = {rTif, aTif, srTif, saTif, cTif};
    float mainNoData[5] = {-2.e9f, -2.e9f, -2.e9f, -2.e9f, -2.e9f};
    makeTiffVRT(cullPar->outFileVRT, mainBands, 5, mainNoData, cullPar->metaData);
    const char *mtBands[1] = {tTif};
    float mtNoData[1] = {0.f};
    makeTiffVRT(cullPar->outFileVRTMT, mtBands, 1, mtNoData, cullPar->metaDataMT);
    free(rTif); free(aTif); free(srTif); free(saTif); free(cTif); free(tTif);
}

void writeCullVrt(CullParams *cullPar, int32_t *byteOrder)
{
    char MTVrt[2048], *tmp, *byteSwapOption;
    // Band stuff for MT files
    GDALDataType dataTypesM[] = {GDT_Byte};
    fprintf(stderr, "writing cull vrt...\n");
    sprintf(MTVrt, "%s.vrt", cullPar->outFileT);
    char *bandFilesM[] = {cullPar->outFileT};
    char *bandNamesM[] = {"MatchType"};
    
    fprintf(stderr, "%s\n", get_value(cullPar->metaData, "ByteOrder"));
    // Set byte order option
    if (strstr(get_value(cullPar->metaData, "ByteOrder"), "MSB"))
    {
        byteSwapOption = "BYTEORDER=MSB";
        *byteOrder = MSB;
    }
    else
    {
        byteSwapOption = "BYTEORDER=LSB";
        *byteOrder = LSB;
    }
    // Band info
    GDALDataType dataTypes[] = {GDT_Float32, GDT_Float32, GDT_Float32, GDT_Float32, GDT_Float32};
    char *bandFiles[] = {cullPar->outFileR, cullPar->outFileA, cullPar->outFileSR, cullPar->outFileSA, cullPar->outFileC};
    char *bandNames[] = {"RangeOffsets", "AzimuthOffsets", "RangeSigma", "AzimuthSigma", "Correlation"};
    // Write VRT
    writeSingleVRT(cullPar->nR, cullPar->nA, cullPar->metaData, cullPar->outFileVRT, bandFiles, bandNames, dataTypes, byteSwapOption, -2.0e9, 5);
    // Write Match VRT
    writeSingleVRT(cullPar->nR, cullPar->nA, cullPar->metaDataMT, MTVrt, bandFilesM, bandNamesM, dataTypesM, NULL, DONOTINCLUDENODATA, 1);
}

static size_t fwriteOptionalBS(void *ptr, size_t nitems, size_t size, FILE *fp, int32_t flags, int32_t byteOrder)
{
    if (byteOrder == LSB)
        return fwrite(ptr, size, nitems, fp);
    else
        return fwriteBS(ptr, nitems, size, fp, flags);
}

void writeCullData(CullParams *cullPar)
{
    FILE *fp;
    int32_t nr, na, byteOrder, i;
    size_t nSamples;
    if (cullPar->tiffFlag) {
        fprintf(stderr, "writing cull data as GeoTIFF...\n");
        writeCullTiffData(cullPar);
        return;
    }
    // Write the vrt
    fprintf(stderr, "writing cull data...\n");
    writeCullVrt(cullPar, &byteOrder);
    nr = cullPar->nR;
    na = cullPar->nA;
    char *files[] = {cullPar->outFileR, cullPar->outFileA, cullPar->outFileSR,
                     cullPar->outFileSA, cullPar->outFileC, cullPar->outFileT};
    void *buffers[] = {cullPar->offRS[0], cullPar->offAS[0], cullPar->sigmaR[0],
                       cullPar->sigmaA[0], cullPar->corr[0], cullPar->type[0]};
    size_t sizes[] = {sizeof(float), sizeof(float), sizeof(float),
                      sizeof(float), sizeof(float), sizeof(char)};
    int32_t flags[] = {FLOAT32FLAG, FLOAT32FLAG, FLOAT32FLAG, FLOAT32FLAG, FLOAT32FLAG, BYTEFLAG};
    nSamples = (cullPar->nA) * cullPar->nR;
    for (i = 0; i < 6; i++)
    {
        fprintf(stderr, "Writing %s\n", files[i]);
        fp = fopen(files[i], "w");
        fwriteOptionalBS(buffers[i], sizes[i], nSamples, fp, flags[i], byteOrder);
        fclose(fp);
    }
 
}
