#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "math.h"
#include <stdlib.h>
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
#include "cullls.h"
static void readArgs(int argc,char *argv[],CullLSParams *cullPar);
static void usage();
/* 
   Global variables definitions
*/
int RangeSize=0;       /* Range size of complex image */
int AzimuthSize=0;   /* Azimuth size of complex image */
int BufferSize=0;     /* Size of nonoverlap region of the buffer */
int BufferLines = 0;         /* # of lines of nonoverlap in buffer */
double RangePixelSize=0; /* Range PixelSize */
double AzimuthPixelSize=0; /* Azimuth PixelSize */
double SLat=-91.;
static void writeLSCulledOffsets(CullLSParams *cullPar);

void  main(int argc, char *argv[])
{   
	CullLSParams cullPar;

	readArgs(argc,argv,&cullPar);
	fprintf(stderr,"cullIslandThresh %i\n",cullPar.islandThresh);
	/* 
	   Load data
	*/
	loadLSCullData(&cullPar);
	/*
	  cull data
	*/
	fprintf(stderr,"Cull data\n");
	cullPar.nMatch=-1; /* Use negative value so this only gets set on the first pass through */
	cullLSData(&cullPar);
	cullLSData(&cullPar);
	cullLSData(&cullPar); 
	/*
	  Compute stats for data
	*/
	fprintf(stderr,"Cull Stats\n");
	cullLSStats(&cullPar);
	/*
	  Compute stats for data
	*/
	fprintf(stderr,"Cull Smooth\n");
	cullLSSmooth(&cullPar);
	/*
	  cull islands
	*/
	if(cullPar.islandThresh > 0) cullLSIslands(&cullPar);
	/*
	  Output result
	*/
	fprintf(stderr,"Cull write\n");
	writeLSCulledOffsets(&cullPar);
}

/******************************************************************************************
   Write culled offsets 
*******************************************************************************************/ 
static void writeLSCulledOffsets(CullLSParams *cullPar)
{
	uint32 i,j;
	FILE *fp;
	char *file1;
	size_t  sl;
	int32 nx,ny;
	extern int32 nMatch, nAttempt, nTotal;
	int32 nGoodPts;
	double sigmaXAvg,sigmaYAvg;
	sl=1500;
	file1=(char *)malloc(sl);

	nx=cullPar->matches.nx;
	ny=cullPar->matches.ny;
	nGoodPts=0;
	/*
	  Final count
	*/
	sigmaXAvg=0; 	sigmaYAvg=0;
	for(i=0; i < ny; i++) {for(j=0; j < nx ; j++) {
			if(cullPar->XS[i][j] > (NODATA+1) &&  cullPar->YS[i][j] > (NODATA+1) )  {
				sigmaXAvg += cullPar->matches.sigmaX[i][j] * cullPar->matches.sigmaX[i][j] ;
				sigmaYAvg += cullPar->matches.sigmaY[i][j] * cullPar->matches.sigmaY[i][j] ;
				nGoodPts++; 
			}
		}
	}
	sigmaXAvg=sqrt(sigmaXAvg/(double) nGoodPts);
	sigmaYAvg=sqrt(sigmaYAvg/(double) nGoodPts);
	/*	fprintf(stderr,"File root %s %i %i\n",matchP->outputFile, nx,ny);
		for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,matchP->outputFile); file1=strcat(file1,".rho"); fprintf(stderr,"%s\n",file1);
		fp=fopen(file1,"w");  fwriteBS(matches->Rho[0],sizeof(float),(size_t)(nx*ny),fp,FLOAT32FLAG); fclose(fp);*/
	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,cullPar->fitDat.matchFile); file1=strcat(file1,".cull.dx"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w");  fwriteBS(cullPar->XS[0],sizeof(float),(size_t)(nx*ny),fp,FLOAT32FLAG); fclose(fp);

	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,cullPar->fitDat.matchFile); file1=strcat(file1,".cull.dy"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w");  fwriteBS(cullPar->YS[0],sizeof(float),(size_t)(nx*ny),fp,FLOAT32FLAG); fclose(fp);

	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,cullPar->fitDat.matchFile); file1=strcat(file1,".cull.sx"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w");  fwriteBS(cullPar->matches.sigmaX[0],sizeof(float),(size_t)(nx*ny),fp,FLOAT32FLAG); fclose(fp);

	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,cullPar->fitDat.matchFile); file1=strcat(file1,".cull.sy"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w");  fwriteBS(cullPar->matches.sigmaY[0],sizeof(float),(size_t)(nx*ny),fp,FLOAT32FLAG); fclose(fp);

	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,cullPar->fitDat.matchFile); file1=strcat(file1,".cull.mtype"); fprintf(stderr,"%s\n",file1);
	fp=fopen(file1,"w");  fwrite(cullPar->matches.type[0],sizeof(uint8),(size_t)(nx*ny),fp); fclose(fp);

	for(i=0; i<sl; i++) file1[i]='\0';  file1=strcpy(file1,cullPar->fitDat.matchFile); file1=strcat(file1,".cull.dat"); fprintf(stderr,"%s\n",file1);

	fp=fopen(file1,"w");
	fprintf(fp,"fileEarly = %s\n",cullPar->matches.fileEarly);
	fprintf(fp,"fileLate = %s\n",cullPar->matches.fileLate);
	fprintf(fp,"x0 = %11.2lf\n",cullPar->matches.x0);
	fprintf(fp,"y0 = %11.2lf\n",cullPar->matches.y0);
	fprintf(fp,"dx = %10.5lf\n",cullPar->matches.dx);

	fprintf(fp,"dy = %10.5lf\n",cullPar->matches.dy);
	fprintf(fp,"stepY = %u\n",cullPar->matches.stepX);
	fprintf(fp,"stepX = %u\n",cullPar->matches.stepY);
	fprintf(fp,"nx = %9i\n",nx);
	fprintf(fp,"ny = %9i\n",ny);
	fprintf(fp,"slowFlag = %9i\n",cullPar->fitDat.slowFlag);
	fprintf(fp,"EPSG = %9i\n",cullPar->fitDat.proj);
	fprintf(fp,"earlyImageJD = %10lf\n",cullPar->matches.jdEarly);
	fprintf(fp,"lateImageJD = %10lf\n",cullPar->matches.jdLate);
	fprintf(fp,"IntervalBetweenImages = %5i\n",(int)(cullPar->matches.jdLate-cullPar->matches.jdEarly));
	fprintf(fp,"Success_rate_for_attempted_matches(%%) =  %7.2f \n",( (float)cullPar->nMatch/(float)cullPar->nAttempt)*100.);
	fprintf(fp,"Culled_rate_for_attempted_matches(%%) =  %7.2f \n",( (float)nGoodPts/(float)cullPar->nAttempt)*100.);
	fprintf(fp,"Mean_sigmaX = %11.3lf\n",sigmaXAvg);
	fprintf(fp,"Mean_sigmaY = %11.3lf\n",sigmaYAvg);
	fprintf(fp,"& \n");
	fclose(fp);


}


/******************************************************************************************
	 read args
*******************************************************************************************/ 
static void readArgs(int argc,char *argv[],CullLSParams *cullPar)
{
	int filenameArg;
	char *argString;
	char *inBase,*outBase;
	int islandThresh;
	int sx,sy;
	int boxSize,nGood;
	float maxY,maxX;
	int i,n,sLen;

	if( argc < 3 || argc > 17 ) usage();        /* Check number of args */ 
	n = argc - 3;
	sx=3;
	sy=3;
	boxSize=9;
	nGood=17;
	maxY=1.0;
	maxX=.75;
	islandThresh=-1;
	for(i=1; i <= n; i++) {
		argString = strchr(argv[i],'-');  
		if(strstr(argString,"sx") != NULL) {
			sscanf(argv[i+1],"%i",&sx);  
			i++;
		} else if(strstr(argString,"sy") != NULL) {
			sscanf(argv[i+1],"%i",&sy);
			i++;
		} else if(strstr(argString,"boxSize") != NULL) {
			sscanf(argv[i+1],"%i",&boxSize);  
			i++;
		} else if(strstr(argString,"nGood") != NULL) {
			sscanf(argv[i+1],"%i",&nGood);  
			i++;
		} else if(strstr(argString,"maxY") != NULL) {
			sscanf(argv[i+1],"%f",&maxY);  
			i++;
		} else if(strstr(argString,"maxX") != NULL) {
			sscanf(argv[i+1],"%f",&maxX);  
			i++;
		} else if(strstr(argString,"islandThresh") != NULL) {
			sscanf(argv[i+1],"%i",&islandThresh);  
			i++;
		} else usage();  
	}
	cullPar->islandThresh=islandThresh;
	cullPar->sX=sx;
	cullPar->sY=sy;
	cullPar->bR=boxSize;
	cullPar->bA=boxSize;
	cullPar->nGood=nGood;
	cullPar->maxY=maxY;
	cullPar->maxX=maxX;
	cullPar->fitDat.matchFile=argv[argc-1];
	return;
}

 
static void usage()
{ 
	error("cullst -islandThresh islandThresh -maxX maxX -maxY maxY -nGood nGood -boxSize boxSize -sx sx -sy sy inbase\n%s\n\n%s\n",
	      "where",
	      "sx,sy = smoothing window size to apply",
	      "islandThresh = cull isolated islands < islandThresh in diameter",
	      "maxX,maxY = max deviation from local median");
}
