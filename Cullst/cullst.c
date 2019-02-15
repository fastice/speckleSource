#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "cullst.h"
#include "math.h"
#include <stdlib.h>
#include <unistd.h>
#define LARGEINT  -2e9
static char **mallocByteMat(int nA,int nR);
static void readArgs(int argc,char *argv[],CullParams *cullPar);
static void usage();
static void addSubtractSimOffsets(CullParams *cullPar,float mySign);
/* 
   Global variables definitions
*/
int RangeSize=0;       /* Range size of complex image */
int AzimuthSize=0;   /* Azimuth size of complex image */
int BufferSize=0;     /* Size of nonoverlap region of the buffer */
int BufferLines = 0;         /* # of lines of nonoverlap in buffer */
double RangePixelSize=0; /* Range PixelSize */
double AzimuthPixelSize=0; /* Azimuth PixelSize */


void main(int argc, char *argv[])
{   
	CullParams cullPar;
	FILE *fp;
    
	readArgs(argc,argv,&cullPar);
  
	fprintf(stderr,"Infiles: \n %s\n %s\n %s\n %s\n %s\n",cullPar.inFileR,cullPar.inFileA,cullPar.inFileC,cullPar.inFileD,cullPar.inFileT);
	fprintf(stderr,"Outfiles: \n %s\n %s\n %s\n %s\n %s\n", cullPar.outFileR,cullPar.outFileA,cullPar.outFileC,   cullPar.outFileD,cullPar.outFileT);
	fprintf(stderr,"cullIslandThresh %i\n",cullPar.islandThresh);
	/* 
	   Load data
	*/
	loadCullData(&cullPar);
	/*
	   if useSim flag set, then read sim offsets, and subtract the intial guess prior to culling.
	 */
	if(cullPar.useSim == TRUE) {
		loadSimData(&cullPar);
		fprintf(stderr,"removing sim offsets \n");		
		addSubtractSimOffsets(&cullPar,-1.0);
	  }
	/*
	  load mask
	*/
	if( access("offsets.mask",F_OK) != -1 && cullPar.ignoreOffsets==FALSE ) {
		/* read mask file if it is the same size, if not skip, and set flag FALSE */
		cullPar.maskFlag=loadCullMask(&cullPar);
	} else cullPar.maskFlag=FALSE;

	/*
	  cull data
	*/
	fprintf(stderr,"Cull data\n");
	cullSTData(&cullPar);
	cullSTData(&cullPar);
	cullSTData(&cullPar);
	/*
	  Compute stats for data
	*/
	fprintf(stderr,"Cull Stats\n");
	cullStats(&cullPar);
  	/*
	   if useSim flag set, then add them back in after the residual difference was culled.
	 */
	if(cullPar.useSim == TRUE) {
		fprintf(stderr,"adding sim offsets back\n");
		 addSubtractSimOffsets(&cullPar,1.0);
	  }	
	/*
	  Compute stats for data
	*/
	fprintf(stderr,"Cull Smooth\n");
	cullSmooth(&cullPar);
	/*
	  cull islands
	*/
	if(cullPar.islandThresh > 0) cullIslands(&cullPar);

	/*
	  Output result
	*/
	fprintf(stderr,"Cull write\n");
	writeCullData(&cullPar);
} 

static void addSubtractSimOffsets(CullParams *cullPar,float mySign) {
	long int i,j;
	for(i=0; i < cullPar->nR; i++)
		for(j=0; j < cullPar->nA; j++) {
			if( cullPar->offA[j][i] < (-LARGEINT+1) && cullPar->offSimA[j][i] < (-LARGEINT+1)) {
				cullPar->offA[j][i] += mySign * cullPar->offSimA[j][i] ;
				cullPar->offR[j][i] += mySign * cullPar->offSimR[j][i] ;				
			}
		}
}             

static void readArgs(int argc,char *argv[],CullParams *cullPar)
{
	int filenameArg;
	char *argString;
	char *inBase,*outBase,*simBase;
	int islandThresh;
	int sr,sa;
	int boxSize,nGood;
	float maxA,maxR;
	int singleMT;
	int i,n,sLen;

	if( argc < 3 || argc > 21 ) {fprintf(stderr,"to Many Argc %i %i\n",argc,21); usage();  }      /* Check number of args */ 
	n = argc - 3;
	sr=2;
	sa=8;
	boxSize=9;
	nGood=17;
	maxA=1.0;
	maxR=.75;
	singleMT=-1;
	islandThresh=-1;
	cullPar->ignoreOffsets=FALSE;
	cullPar->useSim=FALSE;
       
	for(i=1; i <= n; i++) {
		argString = strchr(argv[i],'-');  
		if(strstr(argString,"sr") != NULL) {
			sscanf(argv[i+1],"%i",&sr);  
			i++;
		} else if(strstr(argString,"sa") != NULL) {
			sscanf(argv[i+1],"%i",&sa);
			i++;
		} else if(strstr(argString,"boxSize") != NULL) {
			sscanf(argv[i+1],"%i",&boxSize);  
			i++;
		} else if(strstr(argString,"nGood") != NULL) {
			sscanf(argv[i+1],"%i",&nGood);  
			i++;
		} else if(strstr(argString,"singleMT") != NULL) {
			sscanf(argv[i+1],"%i",&singleMT);
			if( singleMT < 1 ||  singleMT > 3) {
				fprintf(stderr,"singleMT should be in range 0 to 3\n");
				usage();
			}			
			i++;			
		} else if(strstr(argString,"maxA") != NULL) {
			sscanf(argv[i+1],"%f",&maxA);  
			i++;
		} else if(strstr(argString,"maxR") != NULL) {
			sscanf(argv[i+1],"%f",&maxR);  
			i++;
		} else if(strstr(argString,"islandThresh") != NULL) {
			sscanf(argv[i+1],"%i",&islandThresh);  
			i++;
		} else if(strstr(argString,"useSim") != NULL) {
			cullPar->useSim=TRUE;
			simBase=argv[i+1];
			i++;			
		} else if(strstr(argString,"ignoreOffsets") != NULL) {
			cullPar->ignoreOffsets=TRUE;			
		} else { fprintf(stderr,"%i %s\n",i,argv[i]); usage();  }
	}
	cullPar->islandThresh=islandThresh;
	cullPar->sR=sr;
	cullPar->sA=sa;
	cullPar->bR=boxSize;
	cullPar->bA=boxSize;
	cullPar->nGood=nGood;
	cullPar->maxA=maxA;
	cullPar->maxR=maxR;

	cullPar->singleMT=singleMT;
	
	inBase=argv[argc-2];
	outBase=argv[argc-1];

	fprintf(stderr,"|%s|\n",inBase);
	fprintf(stderr,"|%s|\n",outBase);
	sLen = strlen(outBase);
	cullPar->outFileR= (char *)malloc(sizeof(char)*sLen + 4);
	cullPar->outFileA= (char *)malloc(sizeof(char)*sLen + 4);
	cullPar->outFileC= (char *)malloc(sizeof(char)*sLen + 4);
	cullPar->outFileT= (char *)malloc(sizeof(char)*sLen + 4);
	cullPar->outFileD= (char *)malloc(sizeof(char)*sLen + 5);
	sLen = strlen(inBase);
	cullPar->outFileSA= (char *)malloc(sizeof(char)*sLen + 4);
	cullPar->outFileSR= (char *)malloc(sizeof(char)*sLen + 4);

	cullPar->outFileR[0]='\0';
	cullPar->outFileR[0]='\0';
	cullPar->outFileC[0]='\0';
	cullPar->outFileT[0]='\0';
	cullPar->outFileD[0]='\0';
	cullPar->outFileSA[0]='\0';
	cullPar->outFileSR[0]='\0';


	strcat(cullPar->outFileR,outBase);
	strcat(cullPar->outFileA,outBase);
	strcat(cullPar->outFileC,outBase);
	strcat(cullPar->outFileT,outBase);
	strcat(cullPar->outFileD,outBase);
	strcat(cullPar->outFileSA,inBase);
	strcat(cullPar->outFileSR,inBase);

	strcat(cullPar->outFileR,".dr");
	strcat(cullPar->outFileA,".da");
	strcat(cullPar->outFileC,".cc");
	strcat(cullPar->outFileT,".mt");
	strcat(cullPar->outFileD,".dat"); 
	strcat(cullPar->outFileSR,".sr");
	strcat(cullPar->outFileSA,".sa"); 

	sLen = strlen(inBase);
	cullPar->inFileR= (char *)malloc(sizeof(char)*sLen + 4);
	cullPar->inFileA= (char *)malloc(sizeof(char)*sLen + 4);
	cullPar->inFileC= (char *)malloc(sizeof(char)*sLen + 4);
	cullPar->inFileT= (char *)malloc(sizeof(char)*sLen + 4);
	cullPar->inFileD= (char *)malloc(sizeof(char)*sLen + 5);
	cullPar->inFileR[0]='\0';
	cullPar->inFileA[0]='\0';
	cullPar->inFileC[0]='\0';
	cullPar->inFileT[0]='\0';
	cullPar->inFileD[0]='\0';

	strcat(cullPar->inFileR,inBase);
	strcat(cullPar->inFileA,inBase);
	strcat(cullPar->inFileC,inBase);
	strcat(cullPar->inFileT,inBase);
	strcat(cullPar->inFileD,inBase);

	strcat(cullPar->inFileR,".dr");
	strcat(cullPar->inFileA,".da");
	strcat(cullPar->inFileC,".cc");
	strcat(cullPar->inFileT,".mt");
	strcat(cullPar->inFileD,".dat");
       
	if(cullPar->useSim == TRUE) {
		sLen = strlen(simBase);
		cullPar->inFileSimR= (char *)malloc(sizeof(char)*sLen + 4);
		cullPar->inFileSimA= (char *)malloc(sizeof(char)*sLen + 4);
		cullPar->inFileSimD= (char *)malloc(sizeof(char)*sLen + 5);
		cullPar->inFileSimR[0]='\0';
		cullPar->inFileSimA[0]='\0';
		cullPar->inFileSimD[0]='\0';
		strcat(cullPar->inFileSimR,simBase);
		strcat(cullPar->inFileSimA,simBase);
		strcat(cullPar->inFileSimD,simBase);		
		strcat(cullPar->inFileSimR,".dr");
		strcat(cullPar->inFileSimA,".da");
		strcat(cullPar->inFileSimD,".dat");
		fprintf(stderr,"Using simfile : %s %s %s\n",cullPar->inFileSimR,cullPar->inFileSimA,cullPar->inFileSimD);
	}

	return;
}

 
static void usage()
{ 
	error("cullst -useSim offsets -ignoreOffsets -islandThresh islandThresh -singleMT singleMT -maxR maxR -maxA maxA -nGood nGood -boxSize boxSize -sr sr -sa sa inbase outbase\n%s\n\n%s\n%s\n%s\n%s\n%s\n%s\n",
	      "where",
	      "useSim is a simulated offsets file to subtract off in cull",
	      "sr,sa = smoothing window size to apply (if odd then 111, if even 0.51110.5",
	      "islandThresh = cull isolated islands < islandThresh in diameter",
	      "singleMT = only retain this match type (1,2,or3)",
	      "maxR,maxA = max deviation from local median",	      
	      "ignoreOffsets = do not use information from offset file to cull large matches"	      );
}
