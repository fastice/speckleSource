#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "strack.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"
/*
   Read initial Offset information
*/

void parseInitialOffsets(TrackParams *trackPar)
{
    int32_t lineCount = 0, eod;
    char line[256], *tmp;
    FILE *fp;
    int32_t aDone, rDone;
    fprintf(stderr, "\nREADING INTIAL OFFSETS\n");

    if (trackPar->initialOffsetFile != NULL)
        fp = openInputFile(trackPar->initialOffsetFile);
    else
        error("parseInitialOffsets: missing initalOffsetFile\n");
    aDone = FALSE;
    rDone = FALSE;
    while (aDone != TRUE || rDone != TRUE)
    {
        lineCount = getDataString(fp, lineCount, line, &eod);
        if (eod == TRUE)
            error("parseInitialOffsets: missing polynomial");
        if (strstr(line, "range_offset_polynomial:") != NULL)
        {
            tmp = strstr(line, ":");
            tmp++;
            rDone = TRUE;
            if (sscanf(tmp, "%lf%lf%lf", &(trackPar->rShiftPoly[0]),
                       &(trackPar->rShiftPoly[1]), &(trackPar->rShiftPoly[2])) != 3)
                error("parseInitialOffsets: missing range parameter %s", tmp);
        }
        else if (strstr(line, "azimuth_offset_polynomial:") != NULL)
        {
            tmp = strstr(line, ":");
            tmp++;
            if (sscanf(tmp, "%lf%lf%lf", &(trackPar->aShiftPoly[0]),
                       &(trackPar->aShiftPoly[1]), &(trackPar->aShiftPoly[2])) != 3)
                error("parseInitialOffsets: missing azimuth parameter %s", tmp);
            aDone = TRUE;
        }
    } /* Endwhile */

    fprintf(stderr, "rangePoly %f %f %f\n", trackPar->rShiftPoly[0],
            trackPar->rShiftPoly[1], trackPar->rShiftPoly[2]);
    fprintf(stderr, "azimuthPoly %f %f %f\n", trackPar->aShiftPoly[0],
            trackPar->aShiftPoly[1], trackPar->aShiftPoly[2]);
}
