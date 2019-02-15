#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "strack.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"
/*
   Read baseline information
*/

    void parseBase(TrackParams *trackPar)
{   
    int lineCount=0, eod;
    char line[256];
    FILE *fp; 

    if(trackPar->baseFlag == FALSE) return; 

    fprintf(stderr,"\nREADING BASELINE INFORMATION\n");


    if(trackPar->baseParams != NULL) {
          fp = openInputFile(trackPar->baseParams);
        lineCount=getDataString(fp,lineCount,line,&eod);
/*
    Read baseline params
*/
        if( sscanf(line,"%lf%lf%lf%lf",&(trackPar->Bn),&(trackPar->Bp ),
            &(trackPar->dBn ),&( trackPar->dBp)) != 4)
            error("%s  %i","parseBase -- Missing baseline at line:",lineCount);
        trackPar->dBnQ = 0.0;
        trackPar->dBpQ = 0.0;
    } else if(trackPar->baseFile != NULL) {
        error("parseBase: baseFile not implemented yet");
    } else error("parseBase: missing baseline estimate/parameter file");

    fprintf(stderr,"Bn,Bp,dBn,dBp %f %f %f %f\n",trackPar->Bn,trackPar->Bp,
        trackPar->dBn, trackPar->dBp);

} 

