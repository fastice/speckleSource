
#include "strack.h"
#include <libgen.h>

#define STR_BUFFER_SIZE 1024
#define STR_BUFF(fmt, ...) ({                                    \
    char *__buf = (char *)calloc(STR_BUFFER_SIZE, sizeof(char)); \
    snprintf(__buf, STR_BUFFER_SIZE, fmt, ##__VA_ARGS__);        \
    __buf;                                                       \
})


/* Write a flat row-major buffer to a GeoTIFF with a pixel-coord geotransform.
   No vertical flip: row 0 of the buffer becomes the top row of the tif, matching
   the raw-binary / writeSingleVRT convention (and simInSAR's writeFlatTiff). */
static void writeStrackTiff(const char *filename, void *flatData, int32_t width, int32_t height,
                           GDALDataType dataType, float noDataValue, dictNode *metaData)
{
    double pixGT[6] = {-0.5, 1., 0., -0.5, 0., 1.};
    const char *options[] = {"COMPRESS=DEFLATE", NULL};
    GDALDriverH driver = GDALGetDriverByName("GTiff");
    GDALDatasetH ds = GDALCreate(driver, filename, width, height, 1, dataType, (char **)options);
    if (ds == NULL) {
        fprintf(stderr, "writeStrackTiff: cannot create %s\n", filename);
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

/* Build the shared output metadata dictionary (r0/a0/deltas/geodats/windows).
   Used by both the raw-binary VRT and the GeoTIFF paths. */
static dictNode *buildOutputMeta(TrackParams *trackPar)
{
    char geo2[2048], *tmp, filename[2048];
    dictNode *metaData = NULL;
    insert_node(&metaData, "r0", STR_BUFF("%i", trackPar->rStart * trackPar->scaleFactor));
    insert_node(&metaData, "a0", STR_BUFF("%i", trackPar->aStart * trackPar->scaleFactor));
    insert_node(&metaData, "deltaA", STR_BUFF("%i", trackPar->deltaA * trackPar->scaleFactor));
    insert_node(&metaData, "deltaR", STR_BUFF("%i", trackPar->deltaR * trackPar->scaleFactor));
    insert_node(&metaData, "sigmaStreaks", STR_BUFF("%f", 0.));
    insert_node(&metaData, "sigmaRange", STR_BUFF("%f", 0.));
    insert_node(&metaData, "geo1", STR_BUFF("%s", trackPar->intGeodat));
    filename[0] = '\0';
    strcpy(filename, trackPar->imageFile2);
    geo2[0] = '\0';
    tmp = strcat(geo2, dirname(filename));
    tmp = strcat(geo2, "/");
    tmp = strcat(geo2, basename(trackPar->intGeodat));
    insert_node(&metaData, "geo2", STR_BUFF("%s", geo2));
    insert_node(&metaData, "Image1", STR_BUFF("%s", trackPar->imageFile1));
    insert_node(&metaData, "Image2", STR_BUFF("%s", trackPar->imageFile2));
    insert_node(&metaData, "mask", STR_BUFF("%s", trackPar->maskFile));
    insert_node(&metaData, "wR", STR_BUFF("%i", trackPar->wR));
    insert_node(&metaData, "wA", STR_BUFF("%i", trackPar->wA));
    insert_node(&metaData, "wRa", STR_BUFF("%i", trackPar->wRa));
    insert_node(&metaData, "wAa", STR_BUFF("%i", trackPar->wAa));
    insert_node(&metaData, "scaleFactor", STR_BUFF("%i", trackPar->scaleFactor));
    return metaData;
}

/*
   Write offsets as GeoTIFF bands (Range/Azimuth/Correlation + MatchType),
   wrapped by the same .vrt / .mt.vrt the raw path would write. The band tiffs
   are named .RangeOffsets/.AzimuthOffsets/... so makeTiffVRT's filename-derived
   band Descriptions still satisfy readBothOffsetsStrackVrt's "Range"/"Azimuth"
   checks. The streamed raw band files are removed so the output is tiff-only.
 */
void writeTiffFile(TrackParams *trackPar)
{
    char root[2048], rTif[2048], aTif[2048], cTif[2048], tTif[2048], azdTif[2048], *dot;
    dictNode *metaData = buildOutputMeta(trackPar);
    printDictionary(metaData);
    /* Derive the output root by stripping the trailing ".vrt". */
    strcpy(root, trackPar->vrtFile);
    dot = strrchr(root, '.');
    if (dot != NULL) {
        *dot = '\0';
    }
    snprintf(rTif, sizeof(rTif), "%s.RangeOffsets.tif", root);
    snprintf(aTif, sizeof(aTif), "%s.AzimuthOffsets.tif", root);
    snprintf(cTif, sizeof(cTif), "%s.Correlation.tif", root);
    snprintf(tTif, sizeof(tTif), "%s.MatchType.tif", root);
    /* Azimuth defocus keeps its own name rather than a band role, because it
       is deliberately not wrapped by any vrt -- nothing reads it, it is a
       diagnostic. <root>.azd.tif is simply the tiff form of <root>.azd. */
    snprintf(azdTif, sizeof(azdTif), "%s.azd.tif", root);
    /* No vertical flip: row 0 = azimuth 0, matching the raw convention. */
    writeStrackTiff(rTif, trackPar->offR[0], trackPar->nR, trackPar->nA, GDT_Float32, -2.e9f, metaData);
    writeStrackTiff(aTif, trackPar->offA[0], trackPar->nR, trackPar->nA, GDT_Float32, -2.e9f, metaData);
    writeStrackTiff(cTif, trackPar->corr[0], trackPar->nR, trackPar->nA, GDT_Float32, -2.e9f, metaData);
    writeStrackTiff(tTif, trackPar->type[0], trackPar->nR, trackPar->nA, GDT_Byte, 0.f, metaData);
    if (trackPar->aZDefocus != NULL)
    {
        writeStrackTiff(azdTif, trackPar->aZDefocus[0], trackPar->nR, trackPar->nA, GDT_Float32, -2.e9f, metaData);
    }
    /* Wrap: main VRT = Range/Azimuth/Correlation; .mt.vrt = MatchType. */
    const char *mainBands[3] = {rTif, aTif, cTif};
    float mainNoData[3] = {-2.e9f, -2.e9f, -2.e9f};
    makeTiffVRT(trackPar->vrtFile, mainBands, 3, mainNoData, metaData);
    const char *mtBands[1] = {tTif};
    float mtNoData[1] = {0.f};
    makeTiffVRT(trackPar->MTvrtFile, mtBands, 1, mtNoData, metaData);
    /* Remove the streamed raw band files so only the tiffs + vrt remain. */
    remove(trackPar->outFileR);
    remove(trackPar->outFileA);
    remove(trackPar->outFileC);
    remove(trackPar->outFileT);
    if (trackPar->outFileAzDefocus != NULL)
    {
        remove(trackPar->outFileAzDefocus);
    }
    free_dictionary(metaData);
}

/*
   Write offsets data file
 */
void writeVrtFile(TrackParams *trackPar)
{
    char geo2[2048], *tmp, *byteSwapOption, filename[2048];
    dictNode *metaData = NULL;
    if (trackPar->tiffFlag) {
        writeTiffFile(trackPar);
        return;
    }
    /*
     Define Meta Data
    */
    insert_node(&metaData, "r0", STR_BUFF("%i", trackPar->rStart * trackPar->scaleFactor));
    insert_node(&metaData, "a0", STR_BUFF("%i", trackPar->aStart * trackPar->scaleFactor));
    insert_node(&metaData, "deltaA", STR_BUFF("%i", trackPar->deltaA * trackPar->scaleFactor));
    insert_node(&metaData, "deltaR", STR_BUFF("%i", trackPar->deltaR * trackPar->scaleFactor));
    insert_node(&metaData, "sigmaStreaks", STR_BUFF("%f", 0.));
    insert_node(&metaData, "sigmaRange", STR_BUFF("%f", 0.));
    insert_node(&metaData, "geo1", STR_BUFF("%s", trackPar->intGeodat));
    // Copy filename so dirname does not corrupt
    filename[0] = '\0';
	strcpy(filename, trackPar->imageFile2);
    geo2[0] = '\0';
    tmp = strcat(geo2, dirname(filename));
    tmp = strcat(geo2, "/");
    tmp = strcat(geo2, basename(trackPar->intGeodat));
    fprintf(stderr, "%s %s \n", tmp,  basename(trackPar->intGeodat));
    insert_node(&metaData, "geo2", STR_BUFF("%s", geo2));
    // Optional stuff
    insert_node(&metaData, "Image1", STR_BUFF("%s", trackPar->imageFile1));
    insert_node(&metaData, "Image2", STR_BUFF("%s", trackPar->imageFile2));
    insert_node(&metaData, "mask", STR_BUFF("%s", trackPar->maskFile));
    insert_node(&metaData, "wR", STR_BUFF("%i", trackPar->wR));
    insert_node(&metaData, "wA", STR_BUFF("%i", trackPar->wA));
    insert_node(&metaData, "wRa", STR_BUFF("%i", trackPar->wRa));
    insert_node(&metaData, "wAa", STR_BUFF("%i", trackPar->wAa));
    insert_node(&metaData, "scaleFactor", STR_BUFF("%i", trackPar->scaleFactor));
    printDictionary(metaData);
    //fprintf(stderr,"Image 2 %s\n", trackPar->imageFile2);
    // Band info 
    GDALDataType dataTypesM[] = {GDT_Byte};
    char *bandFilesM[] = {trackPar->outFileT};
    char *bandNamesM[] = {"MatchType"};
    writeSingleVRT(trackPar->nR, trackPar->nA, metaData, trackPar->MTvrtFile, bandFilesM, bandNamesM, dataTypesM, NULL, DONOTINCLUDENODATA, 1);
    // Add the byte order for the fp files
    if (trackPar->byteOrder == LSB) {
        byteSwapOption = "BYTEORDER=LSB";
        insert_node(&metaData, "ByteOrder", "LSB");
    }
    else {
        byteSwapOption = "BYTEORDER=MSB"; 
        insert_node(&metaData, "ByteOrder", "MSB");
    }
    // Band info 
    GDALDataType dataTypes[] = {GDT_Float32, GDT_Float32, GDT_Float32};
    char *bandFiles[] = {trackPar->outFileR, trackPar->outFileA, trackPar->outFileC, trackPar->outFileT};
    char *bandNames[] = {"RangeOffsets", "AzimuthOffsets", "Correlation", "MatchType"};
    // Write the VRT
    writeSingleVRT(trackPar->nR, trackPar->nA, metaData, trackPar->vrtFile, bandFiles, bandNames, dataTypes, byteSwapOption, -2.e9, 3);
}
