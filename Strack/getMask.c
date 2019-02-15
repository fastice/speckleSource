#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "strack.h"
#include "math.h"
/*
  Malloc space, read params for, and input complex interferogram
*/

void getMask(TrackParams *trackPar)
{   
	inputImageStructure inputImage;
	char *buf, *file,*datFile,buf1[2048];
	FILE *fp;
	int i,j;
	int r0,a0,nr,na;
	float deltaR,deltaA;
	int lineCount, eod;
	char line[1024];
	
	if(trackPar->maskGeodat==NULL) error("getMask: Missing geodat filename");
	inputImage.stateFlag=TRUE;
	if(trackPar->maskType == GEODATMASK  ) {
		parseInputFile(trackPar->maskGeodat, &inputImage) ;
		trackPar->maskDat.r0=0;
		trackPar->maskDat.a0=0;
		trackPar->maskDat.nr = inputImage.rangeSize; 
		trackPar->maskDat.na = inputImage.azimuthSize; 
		trackPar->maskDat.nrl = inputImage.nRangeLooks; 
		trackPar->maskDat.nal= inputImage.nAzimuthLooks;
	} else if(trackPar->maskType == OFFSETMASK ) {
		file=trackPar->maskGeodat;
 		fp = openInputFile(file);
		lineCount=0;
		lineCount=getDataString(fp,lineCount,line,&eod);
		if( sscanf(line,"%i%i%i%i%f%f",&r0,&a0,&nr,&na,&deltaR,&deltaA) != 6)
			error("%s  %i of %s", "readOffsets -- Missing image parameters at line:",  lineCount,datFile); 
		fclose(fp);
		trackPar->maskDat.r0=r0;
		trackPar->maskDat.a0=a0;
		trackPar->maskDat.nr = nr;
		trackPar->maskDat.na = na;
		trackPar->maskDat.nrl = (int)deltaR;	
		trackPar->maskDat.nal= (int)deltaA;
	} else error("Invalid mask type");
	fprintf(stderr,"\nmask params r0,a0,nr,na,deltaR,deltaA %i %i %i %i %i %i\n",
	trackPar->maskDat.r0,trackPar->maskDat.a0,trackPar->maskDat.nr,trackPar->maskDat.na,
	trackPar->maskDat.nrl,trackPar->maskDat.nal);
	
	trackPar->maskDat.mask = (char **)  malloc(trackPar->maskDat.na * sizeof(char *));
	buf = (char *)  malloc(trackPar->maskDat.na * trackPar->maskDat.nr *sizeof(char));

	for(i=0; i < trackPar->maskDat.na; i++) 
		trackPar->maskDat.mask[i]=&(buf[i*trackPar->maskDat.nr]);

	fp = fopen(trackPar->maskFile,"r");
	fprintf(stderr,"Reading file %s\n",trackPar->maskFile);
	freadBS(buf,sizeof(char),trackPar->maskDat.nr*trackPar->maskDat.na,fp,BYTEFLAG);
   
	for(i=0; i < trackPar->maskDat.na; i+=200) {
		for(j=0; j < trackPar->maskDat.nr; j+=100) fprintf(stderr,"%1i",trackPar->maskDat.mask[i][j]);
		fprintf(stderr,"   %10i\n", i*trackPar->maskDat.nal  + trackPar->maskDat.a0);
	}

	fclose(fp);
} 







