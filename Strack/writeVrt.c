
#include "strack.h"
#include <libgen.h>

#define STR_BUFFER_SIZE 1024
#define STR_BUFF(fmt, ...) ({                                    \
    char *__buf = (char *)calloc(STR_BUFFER_SIZE, sizeof(char)); \
    snprintf(__buf, STR_BUFFER_SIZE, fmt, ##__VA_ARGS__);        \
    __buf;                                                       \
})


/*
   Write offsets data file
 */
void writeVrtFile(TrackParams *trackPar)
{
    char geo2[2048], *tmp, *byteSwapOption, filename[2048];
    dictNode *metaData = NULL;
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
