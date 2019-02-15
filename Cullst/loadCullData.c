#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "cullst.h"
#include "math.h"
#include <stdlib.h>
#include "mosaicSource/common/common.h"

static char **mallocByteMat(int nA,int nR);
static float **mallocFloatMat(int nA,int nR);
static void  loadFile(void **data, char * filename,int size , int flag,int nr,int na );
                         
    
void loadCullData(CullParams *cullPar)
{
	FILE *fp; 
	char line[1024];
	int lineCount, eod,i,j;
	int r0,a0,nr,na;
	int deltaA,deltaR;
	/*
	  Read .dat file
	*/
	fp = openInputFile(cullPar->inFileD);
	lineCount=getDataString(fp,lineCount,line,&eod);
	if( sscanf(line,"%i%i%i%i%i%i",&r0,&a0,&nr,&na,&deltaR,&deltaA) != 6)
		error("%s  %i of %s\n%s",
		      "readOffsets -- Missing image parameters at line:",
		      lineCount,cullPar->inFileD,line); 
	fclose(fp);
	cullPar->nR=nr;              cullPar->nA=na;
	cullPar->rStart=r0;          cullPar->aStart=a0;
	cullPar->deltaA = deltaA;    cullPar->deltaR = deltaR;

	/*
	  Mallocldatam
	*/
	cullPar->offR = mallocFloatMat(na,nr);
	cullPar->offA = mallocFloatMat(na,nr);
	cullPar->offRS = mallocFloatMat(na,nr);
	cullPar->offAS = mallocFloatMat(na,nr);
	cullPar->corr = mallocFloatMat(na,nr);
	cullPar->sigmaR = mallocFloatMat(na,nr);
	cullPar->sigmaA = mallocFloatMat(na,nr);
	cullPar->type = mallocByteMat(na,nr);
	/*
	  Load data 
	*/
	loadFile((void **)cullPar->offR,cullPar->inFileR,sizeof(float),nr,na,FLOAT32FLAG);
	loadFile((void **)cullPar->offA,cullPar->inFileA,sizeof(float),nr,na,FLOAT32FLAG);
	loadFile((void **)cullPar->corr,cullPar->inFileC,sizeof(float),nr,na,FLOAT32FLAG);
	loadFile((void **)cullPar->type,cullPar->inFileT,sizeof(char),nr,na,BYTEFLAG);
	for(i=0; i < nr; i++)
		for(j=0; j < na; j++) if(cullPar->offA[j][i] <(-LARGEINT+1)) cullPar->type[j][i]=0;
}

static void loadFile(void **data, char * filename,int size ,int nr,int na ,int flag) {
	FILE *fp;
	fp = openInputFile(filename);
	freadBS(data[0],size,nr*na,fp,flag);
	fclose(fp);
}

void  loadSimData(CullParams *cullPar)
{
	FILE *fp; 
	char line[1024];
	int lineCount, eod,i,j;
	int r0,a0,nr,na;
	int deltaA,deltaR;

	fp = openInputFile(cullPar->inFileD);
	lineCount=getDataString(fp,lineCount,line,&eod);
	if( sscanf(line,"%i%i%i%i%i%i",&r0,&a0,&nr,&na,&deltaR,&deltaA) != 6)
		error("%s  %i of %s\n%s","readOffsets -- Missing image parameters at line:",lineCount,cullPar->inFileD,line);
	if(na != cullPar->nA || nr != cullPar->nR || r0 != cullPar->rStart || a0 != cullPar->aStart || deltaR != cullPar->deltaR || deltaA != cullPar->deltaA )
		error("UseSim only compatable if sim is same size as offset file");
	fclose(fp);
	fprintf(stderr,"Loading sim offsets\n");
	cullPar->offSimR = mallocFloatMat(na,nr);
	cullPar->offSimA = mallocFloatMat(na,nr);
	/* Read offsets */
	loadFile((void **)cullPar->offSimR,cullPar->inFileSimR,sizeof(float),nr,na,FLOAT32FLAG);
	loadFile((void **)cullPar->offSimA,cullPar->inFileSimA,sizeof(float),nr,na,FLOAT32FLAG);
}


int loadCullMask(CullParams *cullPar)
{
	FILE *fp,*fpD;
	char line[1024];
	int lineCount, eod,i,j;
	int r0,a0,nr,na;
	int deltaA,deltaR;
	
	fprintf(stderr,"Found offsets.mask\n");
	fprintf(stderr,"%d %d\n",cullPar->nA,cullPar->nR);
	/* Read offsets.dat and check the file is the same - mostly to avoid mask register offsets - added 08/23/18 */
	fpD = openInputFile("offsets.dat");
	lineCount=getDataString(fpD,lineCount,line,&eod);
	if( sscanf(line,"%i%i%i%i%i%i",&r0,&a0,&nr,&na,&deltaR,&deltaA) != 6)
		error("%s  %i of %s\n%s", "readOffsets -- Missing image parameters at line:", lineCount,"offsets.dat",line); 
	fclose(fpD);
	
	if( cullPar->nA != na || cullPar->nR != nr) {
		fprintf(stderr,"Ignoring mask because it is not the same size: regoffsets ?\n");
		return(FALSE);
	}
	cullPar->mask = mallocByteMat(cullPar->nA,cullPar->nR);
	loadFile((void **)cullPar->mask,"offsets.mask",sizeof(char),nr,na,BYTEFLAG);
	return TRUE;
 }



static char **mallocByteMat(int nA,int nR)
{
	char *tmp, **tmp1;

	int i;
	tmp = malloc(nR*nA*sizeof(char));
	tmp1 = (char **)malloc(nA*sizeof(char *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
}

                         
static float **mallocFloatMat(int nA,int nR)
{
	float *tmp, **tmp1;
	int i;
	size_t size;
	size=(size_t)nR*(size_t)nA*sizeof(float);
	tmp = malloc(size);
	tmp1 = (float **)malloc(nA*sizeof(float *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
}
