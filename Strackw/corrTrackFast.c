#include <stdlib.h>
#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "strackt.h"
#include "rdfSource/rdfRoutines/SRTMrdf.h"
#include "math.h"
#include <sys/types.h>
#include "time.h"

#define NOVER 10
#define NFAST 16
#define INTPATCHSIZE 8
#define LA 3
#define BAD 0
#define CMATCH 1
#define AMPMATCH 2
#define AMPMATCHLARGE 3
#define NBUFFERLINES 2000
#define OS 2  /* Amplitude oversample factor */
#define OSA 2
static void correlateFast( TrackParams *trackPar,double **tmpS1, double **tmpS2,float **dataS,float **dataR,float *rShift,float *aShift,float *cMax);
static int corrMatch( TrackParams *trackPar,double **sqS, double **sqR,double **tmpS1, double **tmpS2,float *rShift1,float *aShift1,float *cMax);
static int getPatches(int r1, int a1, int r2, int a2, FILE *fp1,FILE *fp2,TrackParams *trackPar);
static int getCorrPatchesFast(int r1,int a1,int r2, int a2, FILE *fp1,FILE *fp2,TrackParams *trackPar,int large);
static void findImage2Pos(int r1, int a1, TrackParams *trackPar, int *r2, int *a2);
static void fftCorrPlans(TrackParams *trackPar);
static void mallocSpace(TrackParams *trackPar);
static char **mallocByteMat(int nA,int nR);
static float **mallocFloatMat(int nx, int ny);
static double **mallocDoubleMat(int nA,int nR);
static fftw_complex **mallocfftw_complexMat(int nx,int ny);
static  strackComplex **mallocErs1ComplexMat(int nA,int nR);
static void overSampleC(TrackParams *trackPar);
static void cmpTrackFast(TrackParams *trackPar, int *iMax, int *jMax, double *cMax);
static void ampMatchEdge(TrackParams *trackPar,int *iMax,int *jMax, double *cMax,int large);
static void findCPeak(float **corrF, fftw_complex **corrNoPad, int wR, int wA, int *iMax,int *jMax,double *cMax);
static void writeDatFile(TrackParams *trackPar);
static void getShifts(double range,double azimuth, Offsets *offsets,double *dr, double *da);
static void readBothOffsetsStrackw( Offsets *offsets);
static int maskValue(TrackParams *trackPar,int r1,int a1);


  
/* ************************ GLOBAL VARIABLES ***************************/
strackComplex *cBuf1,*cBuf2;         /* Input buffers */
fftw_complex *cfBuf1,*cfBuf2;         /* Input buffers */
fftw_complex **fftF1,**fftF2;      /* FFt's */
fftw_complex **patch1in, **patch2in; /* ... */
fftw_complex  **psNoPad;
fftw_complex **cFast,**cFastOver;
fftw_complex **cNoPad;
fftw_complex **psFast,**psFastOver;
fftw_complex **f1, **f2;
float **c, **cFastOverMag,**cNoPadMag;
time_t startTime,lastTime;	
/* Amplitude stuff */
fftw_complex **img1, **img2;       /* Detected images for amplitude match*/
fftw_complex **psAmpNoPad;
fftw_complex **fftFa1,**fftFa2;      /* FFt's */
fftw_complex **caNoPad;
float **caNoPadMag;
/* Double Amplitude stuff */
fftw_complex **img1L, **img2L;  /* Detected images for amplitude match*/
fftw_complex **psAmpNoPadL;
fftw_complex **fftFa1L,**fftFa2L;      /* FFt's */
fftw_complex **caNoPadL;
float **caNoPadMagL;
/* Phase corrections */
double *thetaD;
double *sinThetaD;
double *cosThetaD;
double *overRange;
fftw_complex *bnCorrect;
fftw_complex **intPatch;
fftw_complex **fftIntPatch;
fftw_complex **intPatchOver;
fftw_complex **fftIntPatchOver;
/* FFT plans */
fftwnd_plan aForward;
fftwnd_plan aReverseNoPad;
StrackBuf imageBuf1;
StrackBuf imageBuf2;     
/* New stuff for oversampling ^^^ */
fftw_complex **img1in, **img2in;/* ^^^ Detected images for amplitude match*/
fftwnd_plan aForwardIn;
fftw_complex **fftFa1os, **fftFa1Los,**fftFa2os, **fftFa2Los; /*  ^^^ */
fftw_complex **fftF1os,**fftF2os; /* ... */
int pccount=0;
fftw_complex **patch1, **patch2;
double **meanS;  
double **sigmaS;
double **corrResult;
float **dataS;
float **dataR;

/* ********************************************************************/

/*
  Main routine for amplitude matching with an edge pad. 
*/

void corrTrackFast(TrackParams *trackPar)
{
	extern float  **cFastOverMag,**cNoPadMag;
	extern fftw_complex **psNoPad;
	extern int pccount;
	extern double **meanS;  
	extern double **sigmaS;
	extern double **corrResult;
	double **tmpS1,**tmpS2;
	extern float **dataS,**dataR;
	double **sqS,**sqR;
	FILE *fp1,*fp2;
	FILE *fpR, *fpA, *fpC, *fpT; 
	double tmpx1,tmpy1;
	int i,j,nGood;
	int r1,a1,r2,a2;
	int r1a,a1a,r2a,a2a;
	int iShift,jShift;
	int wR2, wA2;
	float cMax;
	int WA2, WR2;
	double wRover2exp, wAover2exp;
	extern time_t startTime,lastTime;	
	float rShift,aShift,rShift1,aShift1;
	int jMax,iMax;
	int nTot;
	int good,type;
	double cAvg;
	int prevShift;
	char *tBuf;
	int maskVal;  /* Flag to determine whether to match */
	int nMask;
	fprintf(stderr,"\nSPECKLE TRACKING\n");
	/* 
	   Internally keep edge pad bigger to acount for over sampling around the peak.
	   Program operates though edge pad is size specified (e.g., peak must be within
	   original edge pad limits.
	*/
	trackPar->edgePadR += NFAST/4;
	trackPar->edgePadA += NFAST/4;
	/*
	  Malloc space
	*/
	mallocSpace(trackPar);
	wR2=(trackPar->wRa-2*trackPar->edgePadR)*OS;
	wA2=(trackPar->wAa-2*trackPar->edgePadA)*OS;
	fprintf(stderr,"wR2,wA2 %i %i\n",wR2/2,wA2/2);
	meanS = mallocDoubleMat(2*trackPar->edgePadA*OS+1,2*trackPar->edgePadR*OS+1); 
	sigmaS = mallocDoubleMat(2*trackPar->edgePadA*OS+1,2*trackPar->edgePadR*OS+1); 
	corrResult = mallocDoubleMat(2*trackPar->edgePadA*OS+1,2*trackPar->edgePadR*OS+1); 
	dataS=mallocFloatMat( trackPar->wAa*OS, trackPar->wRa*OS);
	for(i=0; i < trackPar->wRa*OS; i++) for(j=0; j < trackPar->wAa*OS; j++) dataS[j][i]=0;
	dataR=mallocFloatMat( wA2, wR2);  
	for(i=0; i < wR2; i++) for(j=0; j < wA2; j++) dataR[j][i]=0;
	sqS=mallocDoubleMat( trackPar->wAa*OS,trackPar->wRa*OS);
	tmpS1=mallocDoubleMat(trackPar->wAa*OS,trackPar->wRa*OS);
	tmpS2=mallocDoubleMat(trackPar->wAa*OS,trackPar->wRa*OS);
	sqR=mallocDoubleMat(wA2, wR2);

	fpR = fopen(trackPar->outFileR,"w");
	fpA = fopen(trackPar->outFileA,"w");
	fpC = fopen(trackPar->outFileC,"w");
	fpT = fopen(trackPar->outFileT,"w"); 

	fprintf(stderr,"Plans\n");
	fftCorrPlans(trackPar);
	/*
	  Open input files
	*/
#ifdef HPLARGEFILES    
	fp1 = fopen(trackPar->imageFile1,"r");
	if(fp1 == NULL) error("can't open %s\n",trackPar->imageFile1);
	fp2 = fopen(trackPar->imageFile2,"r");
	if(fp2 == NULL) error("can't open %s\n",trackPar->imageFile2);
#else
	fp1 = openInputFile(trackPar->imageFile1);
	fp2 = openInputFile(trackPar->imageFile2);
#endif
	/*
	  Initializations
	*/ 
	wRover2exp= trackPar->wR/2. * NOVER*OSA;
	wAover2exp= trackPar->wA/2. * NOVER*OSA;
	trackPar->nFail=0.0;
	trackPar->nComplex=0.0;
	trackPar->nAmp=0;
	trackPar->nAmpL=0;
	/*
	  Initial offsets initialization 
	*/
	if(trackPar->polyShift==FALSE) {
		tBuf=(char *)malloc(strlen(trackPar->initOffsetsFile)+6);
		fprintf(stderr,"%i\n",strlen(trackPar->initOffsetsFile)+6);
		tBuf[0]='\0';
		tBuf= strcat(tBuf,trackPar->initOffsetsFile);   
		tBuf = strcat(tBuf,".da");    
		trackPar->initOffsets.file=tBuf;
		tBuf=(char *)malloc(strlen(trackPar->initOffsetsFile)+6);
		tBuf[0]='\0';
		tBuf= strcat(tBuf,trackPar->initOffsetsFile);   
		tBuf = strcat(tBuf,".dr");    
		trackPar->initOffsets.rFile=tBuf; 
		readBothOffsetsStrackw( &(trackPar->initOffsets));
	}
	/* clear buffers before starting */
	for(i=0; i < trackPar->nA; i++) {
		for(j=0; j < trackPar->nR; j++) {
			trackPar->offR[i][j] = (float)(-LARGEINT);
			trackPar->offA[i][j] =(float)(-LARGEINT);
			trackPar->corr[i][j] = 0;
			trackPar->type[i][j] = BAD;
		}
	}
	/*
	  Loop to do matching 
	*/
	cAvg=0.0;
	good=FALSE;
	prevShift=-1;
	nMask=0;
	/* Loop on azimuth */

	for(i=0; i < trackPar->nA; i++) {
		lastTime=time(NULL);
		/* Azimuth coord */
		a1= trackPar->aStart + i*trackPar->deltaA - trackPar->wAa/2;
		/* Loop on range */
		fprintf(stderr,"i %i %i %i\n",i,trackPar->nA,trackPar->nR);
		cAvg=0; nGood=0;
		for(j=0; j < trackPar->nR; j++) {
			rShift=0.0; aShift=0.0;
			r1= trackPar->rStart + j*trackPar->deltaR - trackPar->wRa/2;
			/* Find position for second image */
			findImage2Pos(r1,a1,trackPar,&r2,&a2);
			rShift+=(r2-r1); aShift += (a2-a1);
			cMax=0.;
			/* returns 1 if no mask available, 0 or 1 with mask */
			/*maskVal=maskValue(trackPar,(r1+trackPar->wRa/2)*trackPar->scaleFactor,(a1+trackPar->wAa/2)*trackPar->scaleFactor);*/
			/* Fixed 10/29/18 - corrections to r1/a1 applied in maskValue */
			maskVal=maskValue(trackPar,r1*trackPar->scaleFactor,a1*trackPar->scaleFactor);
			if(maskVal != 0 && maskVal != 1) fprintf(stderr,"maskVal %i\n",maskVal);
			nMask+=maskVal;
			if(maskVal == 1) {
				/* check bounds of start postion for each patch (lower left corner) */
				if(r1 > 0 && r1 < trackPar->imageP1.nSlpR &&  a1 > 0 && a1 < trackPar->imageP1.nSlpA &&
				   r2 > 0 && r2 < trackPar->imageP2.nSlpR && a2 > 0 && a2 < trackPar->imageP2.nSlpA ) {
					/* Read patches for amplitude match */
					getCorrPatchesFast(r1,a1,r2,a2,fp1,fp2,trackPar,FALSE);
					/* Do the  match */
					corrMatch(trackPar, sqS,sqR,tmpS1,tmpS2,&rShift1,&aShift1,&cMax);
				} else cMax=0.0; /* end if r1>0...*/
				/* 
				   Accept or reject match 
				*/
				if(cMax > .00) { 
					cAvg += cMax; nGood++;
					good=TRUE; trackPar->nAmp++; type=AMPMATCH;
					rShift -= rShift1; aShift -= aShift1;
				} else { rShift=-LARGEINT; aShift=-LARGEINT; trackPar->nFail++;} 
				/* Save values */
				trackPar->offR[i][j] = rShift;
				trackPar->offA[i][j] = aShift;
				if(rShift > -100000) {trackPar->offR[i][j] *=(float)trackPar->scaleFactor;  trackPar->offA[i][j] *=(float)trackPar->scaleFactor; } 
				trackPar->corr[i][j] = cMax;   /*cMax;*/
				trackPar->type[i][j] = type;   /*cMax;*/
			} else {
				trackPar->offR[i][j] = -LARGEINT;  trackPar->offA[i][j] = -LARGEINT;
				trackPar->corr[i][j] = 0;     trackPar->type[i][j] = 0;   
			} /* end if maskVal...*/
		} /* End for j=0... */

		fwriteBS(trackPar->offR[i],sizeof(float),trackPar->nR,fpR,FLOAT32FLAG);
		fwriteBS(trackPar->offA[i],sizeof(float),trackPar->nR,fpA,FLOAT32FLAG);
		fwriteBS(trackPar->corr[i],sizeof(float),trackPar->nR,fpC,FLOAT32FLAG);
		fwriteBS(trackPar->type[i],sizeof(char),trackPar->nR,fpT, BYTEFLAG);
		nTot=(i+1)*trackPar->nR;
		if(nGood > 0) cAvg=cAvg/(double)(nGood); else cAvg=0;
		fprintf(stderr,
			"\r%6i nTot %i, nMatch %8i %4.1f nComplex %8i %4.1f "
			"nAmp %8i %4.1f (%7i-%7i) nFail %7i %4.1f cAvg(line) %4.2f %i -- %5i --",a1,
			nTot,trackPar->nComplex + trackPar->nAmp +trackPar->nAmpL,100.*(double)
			(trackPar->nComplex+trackPar->nAmp+trackPar->nAmpL)/(double)nTot,
			trackPar->nComplex,100.*(double)(trackPar->nComplex)/(double)(nTot),
			trackPar->nAmp+trackPar->nAmpL,    
			100.*(double)(trackPar->nAmp+trackPar->nAmpL)/(double)(nTot),
			trackPar->nAmp,trackPar->nAmpL,
			trackPar->nFail,    100.*(double)(trackPar->nFail)    /(double)(nTot),
			cAvg,nMask,(int) (time(NULL)-lastTime));
	} /* End for i=0... */
	fclose(fpR);
	fclose(fpA);
	fclose(fpC);
	fclose(fpT);
	writeDatFile(trackPar);
}


 

/*
  Read patches
  For simplicity, this program is mostly a repeat of getAmpPatches. 
  It reads the full size patch for image 2, so that subsequent programs
  will just read the relevant patch when needed. 
*/
static int getCorrPatchesFast(int r1,int a1,int r2, int a2, FILE *fp1,FILE *fp2, TrackParams *trackPar,int large)
{
	/*    long long offset1, offset2;*/
	off_t offset1, offset2;
	int i,j,i1,j1,i2,j2,ja,j2a;
	int m,n,k,l;
	int la1,la2,lr1,lr2; 
	int dOff1,dOff2;
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;     
	extern fftw_complex **img1,**img2;
	extern fftw_complex **img1in, **img2in;/* ^^^ Detected images for amplitude match*/
	extern fftwnd_plan aForwardIn; /* ^^^ */
	extern fftwnd_plan aForwardInL; /* ^^^ */
	extern fftw_complex **fftFa1os, **fftFa1Los,**fftFa2os, **fftFa2Los; /*  ^^^ */
	fftw_complex **im1,**im2,**im1in,**im2in;
	extern strackComplex *cBuf1,*cBuf2; 
	extern fftw_complex *cfBuf1,*cfBuf2; 
	extern fftw_complex **fftFa1,**fftFa2;      /* FFt's */
	fftw_complex **f1,**f2,**f1a,**f2a;
	int azShift,rangeShift;
	strackComplex **ers1Buf1,**ers1Buf2;
	fftw_complex **fftwBuf1,**fftwBuf2;
	int edge;
	int wR2,wA2;
	size_t sSize, wRa,wAa,s1,s2,a1a,a2a;
	int ia,i2a;
	fftw_complex zero;
	float scale;
	void *tmp;
	int tmp1;
	FILE *fpJunk,*fpDebug;
	zero.re=0.0; zero.im=0.0;
	if(a1 >= trackPar->imageP1.nSlpA || a2 >= trackPar->imageP2.nSlpA) {
		error("getpatches");
		return(FALSE);
	} 
	/* set points to buffer for image 1 and image 2*/
	if(trackPar->floatFlag ==TRUE) {
		sSize=sizeof(fftw_complex);
		fftwBuf1 = (fftw_complex **)imageBuf1.buf; fftwBuf2 = (fftw_complex **)imageBuf2.buf;
	}else { 
		sSize=sizeof(strackComplex);
		ers1Buf1 = (strackComplex **)imageBuf1.buf; ers1Buf2 = (strackComplex **)imageBuf2.buf;
	}
	wAa=trackPar->wAa;   wRa=trackPar->wRa;    
	im1=img1; im2=img2; 
	im1in=img1in; im2in=img2in;
	f1=fftFa1os;    f2=fftFa2os;
	f1a=fftFa1;     f2a=fftFa2;
	/* size for search chip */
	wR2=(trackPar->wRa-2*trackPar->edgePadR)*OS;
	wA2=(trackPar->wAa-2*trackPar->edgePadA)*OS;
	/* Added navg values to scale computation to avoid overflow that was occurring with tsx data */
	scale=1./((float)wRa*(float)wAa*(float)wRa*(float)wAa*(float)(OSA*OSA*OSA*OSA) * (OS*OS*(trackPar->navgA+1) * (trackPar->navgR+1))); /* rough guess at scale to avoid fp overflow */
	/*
	  Load buffers if needed (only required if patch is outside buffer)
	*/
	if(a1 < imageBuf1.firstRow || (a1+wAa) > imageBuf1.lastRow) {
		a1a=max(0,a1-trackPar->wAa);
		/*a1a=max(0,a1-trackPar->wRa);****/
		s1=trackPar->imageP1.nSlpR*min(trackPar->imageP1.nSlpA-a1a,NBUFFERLINES);
		offset1 = ((off_t)a1a * trackPar->imageP1.nSlpR)*(off_t)sSize;
		fseeko(fp1,offset1,SEEK_SET);
		if(trackPar->floatFlag ==TRUE) freadBS((void *)imageBuf1.buf[0],sSize,s1,fp1,FLOAT32FLAG);
		else freadBS((void *)imageBuf1.buf[0],sSize,s1,fp1,INT16FLAG);
		imageBuf1.firstRow=a1a;
		imageBuf1.lastRow = a1a + min(trackPar->imageP1.nSlpA-a1a,NBUFFERLINES)-1;
	}

	if(a2 < imageBuf2.firstRow || (a2+wAa) > imageBuf2.lastRow) {
		a2a=max(0,a2-trackPar->wAa);		
		/* a2a=max(0,a2-trackPar->wRa); ****/
		s2=trackPar->imageP2.nSlpR*min(trackPar->imageP2.nSlpA-a2a,NBUFFERLINES);
		offset2 = ((off_t)a2a * trackPar->imageP2.nSlpR)*(off_t)sSize;
		fseeko(fp2,offset2,SEEK_SET);
		if(trackPar->floatFlag ==TRUE) freadBS((void *)imageBuf2.buf[0],sSize,s2,fp2,FLOAT32FLAG);
		else  freadBS((void *)imageBuf2.buf[0],sSize,s2,fp2,INT16FLAG);
		imageBuf2.firstRow=a2a;
		imageBuf2.lastRow = a2a + min(trackPar->imageP2.nSlpA-a2a,NBUFFERLINES)-1;
	}
	/*
	  Get complex data
	*/
	for(i=0; i < wAa; i++) {
		i1= a1-imageBuf1.firstRow + i;
		i2= a2-imageBuf2.firstRow + i;
		if(i1 > NBUFFERLINES || i2 > NBUFFERLINES) error("Buffer exceeded in getCorrPatchesFast \n");
		for(j=0; j < wRa; j++) {
			j1=r1+j;
			j2=r2+j;
			if(trackPar->floatFlag ==TRUE) {
				im1in[i][j].re = (fftw_real)fftwBuf1[i1][j1].re;
				im2in[i][j].re = (fftw_real)fftwBuf2[i2][j2].re;
				im1in[i][j].im = (fftw_real)fftwBuf1[i1][j1].im;
				im2in[i][j].im = (fftw_real)fftwBuf2[i2][j2].im;
			} else {
				im1in[i][j].re = ((fftw_real)ers1Buf1[i1][j1].r) ;
				im2in[i][j].re = ((fftw_real)ers1Buf2[i2][j2].r) ;
				im1in[i][j].im = ((fftw_real)ers1Buf1[i1][j1].i);
				im2in[i][j].im = ((fftw_real)ers1Buf2[i2][j2].i);
			}
		}
	}
	/* Forward FFT */
	fftwnd_one(aForwardIn,im1in[0],f1[0]);
	fftwnd_one(aForwardIn,im2in[0],f2[0]);
	/*
	  Zero pad fft for over sampling
	*/
	for(i=0; i < wAa*OSA; i++ ) {
		for(j=0; j < wRa*OSA; j++ ) {
			f1a[i][j]=zero;
			f2a[i][j]=zero;
		}
	}
	azShift=0; rangeShift=0;
	for(i=0; i < wAa/2; i++ ) {
		i1=OSA*wAa - wAa/OSA + i;
		i2=(wAa/OSA + i + azShift) % wAa;
		ia=(i+azShift) % wAa;
		j1=wRa*OSA - wRa/OSA;
		j2= wRa/OSA + rangeShift;
		for(j=0; j < wRa/2; j++ ) {
			ja=(j+rangeShift) % wRa;
			j2a = j2 % wRa;
			f1a[i][j]  = f1[ia][ja];   f1a[i1][j] = f1[i2][ja];
			f1a[i1][j1]= f1[i2][j2a];  f1a[i][j1] = f1[ia][j2a];
			f2a[i][j]  = f2[ia][ja];   f2a[i1][j] = f2[i2][ja];
			f2a[i1][j1]= f2[i2][j2a];  f2a[i][j1] = f2[ia][j2a];
			j1++; j2++;
		} 
	}
	/*
	  Inverse transform
	*/
	fftwnd_one(aReverseNoPad,f1a[0],im1[0]);
	fftwnd_one(aReverseNoPad,f2a[0],im2[0]);
	/*
	  Detect Data , Note the image patches are double bookkept as im1/im2 and dataR, dataS
	  based on the way the program was kluged together
	*/
	/* 
	   Smaller patch for second image 
	*/
	for(i=0; i < wA2; i++) {
		i1=trackPar->edgePadA*OS +i;
		for(j=0; j < wR2; j++) {
			j1=trackPar->edgePadR*OS +j;
			/* Smoothing */			
			if(trackPar->navgA > 1 || trackPar->navgR > 1) {
				dataR[i][j]=0;
				for( m=-trackPar->navgA; m <= trackPar->navgA; m++) {
					k=min(max(i1+m,0),trackPar->edgePadR*OS+wA2+trackPar->navgA); /* Changed 11/12/08 data outside window */
					for( n=-trackPar->navgR; n <= trackPar->navgR; n++) {
						/* fixed 11/12/08, changed from m to n */
						l=min(max(j1+n,0),trackPar->edgePadR*OS+wR2+trackPar->navgR);
						dataR[i][j]+=(im2[k][l].re * im2[k][l].re + im2[k][l].im * im2[k][l].im)*scale;
					}
				}
			} else dataR[i][j]=(im2[i1][j1].re * im2[i1][j1].re + im2[i1][j1].im * im2[i1][j1].im)*scale;
		}
	}
	/*	FILE *fpDebug; 	fpDebug= fopen("test1.debug","w");	fwriteBS(dataR[0],sizeof(float),wR2*wA2,fpDebug,FLOAT32FLAG);        exit(-1);	*/
	/*
	  Full patch
	*/ 
	for(i=0; i < wAa*OS; i++) {
		for(j=0; j < wRa*OS; j++) {
			/* Smoothing */
			if(trackPar->navgA > 1 || trackPar->navgR > 1) {
				dataS[i][j]=0;
				for( m=-trackPar->navgA; m <= trackPar->navgA; m++) {
					k=min(max(i+m,0),wAa*OS-1);
					for( n=-trackPar->navgR; n <= trackPar->navgR; n++) {
						/* fixed 11/12/09, changed from m to n */
						l=min(max(j+n,0),wRa*OS-1);
						dataS[i][j] += (im1[k][l].re * im1[k][l].re + im1[k][l].im * im1[k][l].im)*scale;
					}
				}
			} else {
				dataS[i][j]=(im1[i][j].re * im1[i][j].re + im1[i][j].im * im1[i][j].im)*scale;
			}
		} /* end j */
	} /* end i */
	/*
	  copy  dataS,datR back to im1 and im2 - added sqrt 10/23/18
	*/
	/*fprintf(stderr,"===M\n");	*/
	for(i=0; i < wAa*OS; i++) {
		for(j=0; j < wRa*OS; j++) {
			dataS[i][j]=(float)sqrt((double)dataS[i][j]);
			im1[i][j].re=dataS[i][j] ;
			im1[i][j].im=0.0;
			im2[i][j].im=0.0;
			im2[i][j].re=0.0;
		}
	}
	/*fprintf(stderr,"===N\n");*/
	for(i=0; i < wA2; i++) {
		i1=trackPar->edgePadA*OS +i;
		for(j=0; j < wR2; j++) {
			j1=trackPar->edgePadR*OS +j;
			dataR[i][j]=(float)sqrt((double)dataR[i][j]);
			im2[i1][j1].re=dataR[i][j];
		}
	}
	/*
	fpDebug= fopen("test1.debug","w");	fwriteBS(im1[0],sizeof(float),2*wAa*wRa,fpDebug,FLOAT32FLAG);   
	fpDebug= fopen("test2.debug","w");	fwriteBS(im2[0],sizeof(float),2*wAa*wRa,fpDebug,FLOAT32FLAG);    exit(-1);	
	*/
	return(TRUE);
}



/*************************************************************************
   Find inital shift
   ************************************************************************* */
static void findImage2Pos(int r1, int a1, TrackParams *trackPar, 
			  int *r2, int *a2)
{
	double aShift,rShift;
	double sc;
	sc=(double)trackPar->scaleFactor;
	if(trackPar->polyShift==TRUE) {
		rShift = trackPar->rShiftPoly[0] + (double)(r1*sc)* trackPar->rShiftPoly[1] +
			(double)(a1*sc)* trackPar->rShiftPoly[2];
		aShift = trackPar->aShiftPoly[0] + (double)(r1*sc)* trackPar->aShiftPoly[1] +
			(double)(a1*sc)* trackPar->aShiftPoly[2];
	} else {
		getShifts(r1,a1 ,&(trackPar->initOffsets),&rShift,&aShift);
		/* If offset not available use poly - added 6/2/04 */
		if(rShift < (-LARGEINT+1)) {
			rShift = trackPar->rShiftPoly[0] + (double)(r1*sc)* trackPar->rShiftPoly[1] +
				(double)(a1*sc)* trackPar->rShiftPoly[2];
			aShift = trackPar->aShiftPoly[0] + (double)(r1*sc)* trackPar->aShiftPoly[1] +
				(double)(a1*sc)* trackPar->aShiftPoly[2];
		}
	}
	/* Added 8/21 for tsx scaling */
	rShift /=sc; aShift/=sc; 
	if(rShift >= 0) *r2 = r1 + (int)(rShift + 0.5);
	else *r2 = r1 + (int)(rShift - 0.5);
	if(aShift >= 0) *a2 = a1 + (int)(aShift + 0.5);
	else *a2 = a1 + (int)(aShift - 0.5);
}


/***************************************************************************
   Get shifts from previous estimate
****************************************************************************/
static void getShifts(double range,double azimuth,Offsets *offsets,double *dr, double *da)
{
	float t,u;
	float p1,p2,p3,p4;
	float **rimage, **aimage;
	float minvalue;
	float zeroOffset;
	double bn,bp, bSq;
	double deltaZ;
	float alonTrack;
	double imageLength, normAzimuth;
	int i,j;

	minvalue = (float) -LARGEINT * 0.9;
	/*
	  If out of input image, return 0
	*/
	range=(range-offsets->rO)/offsets->deltaR;
	azimuth=(azimuth-offsets->aO)/offsets->deltaA;
	/* Out of bounds */   
	if(range < 0.0 || azimuth < 0.0 ||
	   (int)range >= (offsets->nr) ||
	   (int)azimuth >= (offsets->na) ) {
		*dr=-LARGEINT;  *da=-LARGEINT;
		return;
	}
	rimage = (float **)offsets->dr;  aimage = (float **)offsets->da;
	j = (int)range;    i = (int)azimuth;
	/*
	  Handle border pixels
	*/
	if(j ==  (offsets->nr-1) ||
	   i == (int) (offsets->na-1) ) {
		*dr =  rimage[i][j];  *da =  aimage[i][j];
		if(*dr < minvalue || *da < minvalue) {*dr=-LARGEINT;  *da=-LARGEINT;}
		return;
	}
	t = (float)(range - (double)j);
	u = (float)(azimuth - (double)i);
	/* Interp r */
	p1 = rimage[i][j];        p2 = rimage[i][j+1];
	p3 = rimage[i+1][j+1];    p4 = rimage[i+1][j];
	if(p1 <= minvalue || p2 <= minvalue || p3 <= minvalue || p4 <= minvalue){
		*dr=max(max(max(p1,p2),p3),p4);
		if(*dr < minvalue) {*dr=-LARGEINT; *da=-LARGEINT; return;}
         
	} else {
		*dr = (float)((1.0 - t)*(1.0 - u) * p1 + t * (1.0 - u) * p2 +
			      t * u * p3 +                (1.0 - t) * u* p4);
	}
	/* Interp a */
	p1 = aimage[i][j];        p2 = aimage[i][j+1];
	p3 = aimage[i+1][j+1];    p4 = aimage[i+1][j];
	if(p1 <= minvalue || p2 <= minvalue || p3 <= minvalue || p4 <= minvalue){
		*da=max(max(max(p1,p2),p3),p4);
		if(*da < minvalue) {*dr=-LARGEINT; *da=-LARGEINT; return;}
	} else {
		*da = (float)((1.0 - t)*(1.0 - u) * p1 + t * (1.0 - u) * p2 +
			      t * u * p3 +                (1.0 - t) * u* p4);
	}
}




/***************************************************************************
    read offsets if used as a first guess
***************************************************************************/
static void readBothOffsetsStrackw( Offsets *offsets) 
{
	char *datFile, buf[1024];
	char line[1024];
	int lineCount, eod;
	FILE *fp;
	int rO,aO,nr,na;
	float deltaA,deltaR;
	double c1;
	float *fBuf,*fBufR;
	char *file;
	int i;
	/*
	  Read inputfile
	*/
	file=offsets->file;
	datFile=&(buf[0]);
	buf[0]='\0';

	datFile= strcat(datFile,file);   
	datFile = strcat(datFile,".dat");
 
	fp = openInputFile(datFile);
	lineCount=0;
	lineCount=getDataString(fp,lineCount,line,&eod);

	if( sscanf(line,"%i%i%i%i%f%f",&rO,&aO,&nr,&na,&deltaR,&deltaA) != 6)
		error("%s  %i of %s",
		      "readOffsets -- Missing image parameters at line:",
		      lineCount,datFile); 
	fclose(fp);

	offsets->nr=nr;
	offsets->na=na;
	offsets->rO=rO;
	offsets->aO=aO;
	offsets->deltaA = deltaA;
	offsets->deltaR = deltaR;

	offsets->da = (float **)malloc(sizeof(float *)*na);
	offsets->dr = (float **)malloc(sizeof(float *)*na);
	/*
	  Malloc buffers
	*/
	fBuf=(float *)malloc(sizeof(float)*na*nr);
	fBufR=(float *)malloc(sizeof(float)*na*nr);
	for(i=0; i < na; i++) {
		offsets->da[i]=&(fBuf[i*nr]); offsets->dr[i]=&(fBufR[i*nr]);
	}
	/*
	  Read files
	*/
	fprintf(stderr,"--- Reading offset file %s\n\n",file);
	fp = openInputFile(offsets->file);
	for(i=0; i < na; i++) freadBS(offsets->da[i],sizeof(float),nr,fp,FLOAT32FLAG);
	fclose(fp);

	fprintf(stderr,"--- Reading offset file %s\n\n",offsets->rFile);
	fp = openInputFile(offsets->rFile);
	for(i=0; i < na; i++) freadBS(offsets->dr[i],sizeof(float),nr,fp,FLOAT32FLAG);
	fclose(fp);
}



/***************************************************************
    Main driver for amplitude match (with edgepad)
***************************************************************/
static int corrMatch( TrackParams *trackPar,double **sqS, double **sqR,double **tmpS1, double **tmpS2,float *rShift1,float *aShift1,float *cMax)
{
	extern float **dataS, **dataR;
	int countR,countS; /* Count of good pixel and reference data */
	int ls,ss;
	int lr,sr;
	int i,j;
	int wR2,wA2;
	double  noData;
	/*
	  Precompute squares and count pixels
	  Return if not enough valid pixels.
	*/
	wR2=(trackPar->wRa-2*trackPar->edgePadR)*OS;
	wA2=(trackPar->wAa-2*trackPar->edgePadA)*OS;
	countR = 0; countS=0;

	for(lr=0; lr < trackPar->wAa*OS; lr++)
		for(sr=0; sr < trackPar->wRa*OS; sr++)
			sqS[lr][sr] = dataS[lr][sr] * dataS[lr][sr];

	for(lr=0; lr < wA2; lr++)
		for(sr=0; sr < wR2; sr++) {
			sqR[lr][sr] =  dataR[lr][sr] * dataR[lr][sr];
		}

	correlateFast(trackPar,tmpS1,tmpS2,dataS,dataR,rShift1,aShift1,cMax); 

	return TRUE;
}

/***********************************************************
   Main code to execute point correlation
************************************************************/

static void correlateFast( TrackParams *trackPar,double **tmpS1, double **tmpS2,float **dataS, float **dataR,float *rShift,float *aShift,float *cMax)
{
	extern double **meanS;  
	extern double **sigmaS;
	/*    extern float **dataS, **dataR;*/
	extern fftw_complex **caNoPad;
	extern float **caNoPadMag;
	extern fftw_complex **cFast,**cFastOver;
	double maxCorr,corr;   
	extern double **corrResult;
	FILE *fpDebug;
	int l,s;
	int l1,s1,s2,l2;
	int iMax,jMax,iMax1,jMax1;
	double rShift1,aShift1;
	double cM;
	int lr,la;
	int i,j,i1,j1;
	int m,n;
	int wR2, wA2;
	double sum1,sum2;
	double meanR,sigmaR;
	double scaleW;
	/* size of search chip */
	wR2=(trackPar->wRa-2*trackPar->edgePadR)*OS;
	wA2=(trackPar->wAa-2*trackPar->edgePadA)*OS;
	/*
	  Step 1: Compute mean/sigma for second image patch
	*/
	sigmaR=0;  meanR=0;
	for(la=0; la < wA2; la++) 
		for(lr=0; lr < wR2; lr++) {
			sigmaR += dataR[la][lr] * dataR[la][lr];
			meanR += dataR[la][lr];
		}
	meanR /= (double)(wR2*wA2);
	sigmaR=sigmaR/(double)(wR2*wA2);
	sigmaR-= meanR*meanR;
	/* 
	   Step 2:
	   Compute mean and variance for S1
	   first,do running average over columns
	*/
 	for(l=0; l < trackPar->wAa*OS; l++) {
		sum1=0.0; sum2=0.0;  s1 = 0;
		for(s=0; s < trackPar->wRa*OS; s++) {
			sum1 += dataS[l][s]; sum2 += dataS[l][s]*dataS[l][s];
			if(s >= (wR2-1) ) {
				s2=s-(wR2-1);
				tmpS1[l][s1] = sum1; tmpS2[l][s1] = sum2;
				sum1 -= dataS[l][s2]; sum2 -= dataS[l][s2]*dataS[l][s2];
				s1++;
			}
		}
	}
	/*
	  finish by do running average over rows
	*/  
	s1=0;
	for(s=0;  s < (2*trackPar->edgePadR*OS)+1; s++) {
		sum1=0.0; sum2=0.0;  l1 = 0;
		for(l=0; l < trackPar->wAa*OS; l++) {
			sum1 += tmpS1[l][s];  sum2 += tmpS2[l][s];
			if(l >= (wA2-1)) {
				l2=l-(wA2-1);
				meanS[l1][s1] = sum1/(double)(wR2*wA2);
				sigmaS[l1][s1] = sum2/(double)(wR2*wA2) - meanS[l1][s1]*meanS[l1][s1];
				sum1 -= tmpS1[l2][s];  sum2 -= tmpS2[l2][s];
				l1++;
			}
		}
		s1++;
	}
	/*
	  Step 3, do fft convolution of image patchtes
	*/
	ampMatchEdge(trackPar,&iMax,&jMax,&cM,FALSE); 
	/*
	  Step 4: find peak in correlation function
	*/    
	scaleW=1.0/( (double)(OS*OS*OS*OS)*(double)trackPar->wAa*(double)(trackPar->wAa-2*trackPar->edgePadA) * 
		     (double)trackPar->wRa* (trackPar->wRa - 2*trackPar->edgePadR) ) ;
	maxCorr=-1E30;
	/* Loop to search over.Note goes past edge pad to fill buffer for over sampling */
	for(l=-OS*trackPar->edgePadA; l <= OS*trackPar->edgePadA; l++) {
		i=l+OS*trackPar->edgePadA;
		m=l+trackPar->wAa*OS/2;
		for(s=-2*trackPar->edgePadR; s <= 2*trackPar->edgePadR; s++) {
			j=s+OS*trackPar->edgePadR;
			n=s+trackPar->wRa*OS/2;
			corrResult[i][j] = caNoPadMag[m][n] * scaleW;
			corrResult[i][j] = corrResult[i][j] -  (meanR*meanS[i][j]);
			corrResult[i][j] = corrResult[i][j] / sqrt( sigmaR * sigmaS[i][j] );
			/* Only select point if within original specified edge pad */
			if(corrResult[i][j] > maxCorr && 
			   l >= -OS*(trackPar->edgePadA-NFAST/4) &&
			   s >= -OS*(trackPar->edgePadR-NFAST/4) &&
			   l <=  OS*(trackPar->edgePadA-NFAST/4) &&
			   s <=  OS*(trackPar->edgePadR-NFAST/4)    ) {
				maxCorr=corrResult[i][j];
				iMax=i; jMax=j;
				/*fprintf(stderr,"%f %lf %lf m %lf %lf s %lf %lf \n",caNoPadMag[m][n] ,scaleW,caNoPadMag[m][n] * scaleW,meanR,meanS[i][j],sigmaR,sigmaS[i][j]);*/
			}
		}
	}
	/*fprintf(stderr,"%f %f %f %lf %f %f \n",maxCorr,corrResult[iMax][jMax],meanR,sigmaR,meanS[iMax][jMax],sigmaS[iMax][jMax] );*/
	if(maxCorr < 0) return;
	/*
	  Step 5:  load from around peak to oversample
	*/
	for(i=0; i < NFAST; i++) {
		i1 = iMax -NFAST/2 + i; j1 = jMax - NFAST/2;
		for(j=0; j < NFAST; j++) {
			cFast[i][j].re =corrResult[i1][j1];  cFast[i][j].im = 0.0;
			j1++;
		}
	}
	/*
	  Step 6: Oversample correlatin peak
	*/
	overSampleC(trackPar);
	/*
	  Step 7: Find peak and renorm value
	*/
	*cMax=-1.0e30;
	for(i=(NFAST*NOVER)/2-NOVER*2; i <  (NFAST*NOVER)/2 + NOVER*2; i++)
		for(j=(NFAST*NOVER)/2-NOVER*2; j <  (NFAST*NOVER)/2 + NOVER*2; j++) {
			corr = (cFastOver[i][j].re * cFastOver[i][j].re +cFastOver[i][j].im * cFastOver[i][j].im);
			if(corr > *cMax) {*cMax=corr; iMax1=i; jMax1=j;}   
		}
	*cMax=sqrt(*cMax)/(NFAST*NFAST);
	/*
	  Step 8: Compute raw, fractional pixel shift. 
	*/
	iMax= (iMax * NOVER) + iMax1-(NOVER*NFAST)/2;
	jMax= (jMax * NOVER) + jMax1-(NOVER*NFAST)/2;
	*rShift =((float)jMax-trackPar->edgePadR * OS * NOVER) /(NOVER * OS);
	*aShift =((float)iMax-trackPar->edgePadA * OS * NOVER) /(NOVER * OS);
}



/*********************************************************
   FFT's to do cross correlation for matching
**********************************************************/
static void ampMatchEdge(TrackParams *trackPar,int *iMax,int *jMax,double *cMax,int large)
{
	extern fftwnd_plan aForward;
	extern fftwnd_plan aReverseNoPad;
	extern fftw_complex **psAmpNoPad;
	extern float **caNoPadMag;
	extern fftw_complex **caNoPad;
	extern fftw_complex **fftFa1,**fftFa2;
	extern fftw_complex **img1,**img2;
	fftwnd_plan *aFor;
	fftw_complex **im1,**im2;
	extern fftw_complex **psAmpNoPad;
	fftw_complex **psANoPad;  
	fftw_complex **fa1,**fa2;
	int wR,wA;
	double avgAmp1,avgAmp2;
	double p1,p2;
	int i,j,i1,j1;

	fa1=fftFa1; fa2=fftFa2;
	im1=img1; im2=img2;
	psANoPad = psAmpNoPad;
	wA = trackPar->wAa; wR = trackPar->wRa;
	aFor = &aForward;
	psANoPad = psAmpNoPad;
	/*
	  Step 1: FFT images  power spectrum
	*/
	fftwnd_one(*aFor,im1[0],fa1[0]);
	fftwnd_one(*aFor,im2[0],fa2[0]);    
	/*
	  Step 2: Compute zero padded power spectrum
	*/
	i1=0;
	for(i=0; i < wA*OS; i++ ) {
		for(j=0; j < wR*OS; j++ ) {
			psANoPad[i][j].re = fa1[i1][j].re * fa2[i1][j].re + fa1[i1][j].im * fa2[i1][j].im; 
			psANoPad[i][j].im = fa1[i1][j].im * fa2[i1][j].re - fa1[i1][j].re * fa2[i1][j].im;
		} /* End for j */
		i1++;
	} /* End for i */
	/*
	  Step 3: Inverse transform to get convolved image patches
	*/
	fftwnd_one(aReverseNoPad,psAmpNoPad[0],caNoPad[0]);
	/* 
	   Step 4: Unscramble result
	*/ 
	for(i=0; i < wA*OS/2; i++ ) {
		i1=(wA*OS)/2 + i;
		j1=(wR*OS)/2;
		for(j=0; j < wR*OS/2; j++ ) {
			caNoPadMag[i][j]  	= sqrt(caNoPad[i1][j1].re   * caNoPad[i1][j1].re + caNoPad[i1][j1].im * caNoPad[i1][j1].im);
			caNoPadMag[i1][j1] 	= sqrt(caNoPad[i][j].re	 * caNoPad[i][j].re 	+ caNoPad[i][j].im 	 * caNoPad[i][j].im);
			caNoPadMag[i1][j]  	= sqrt(caNoPad[i][j1].re 	* caNoPad[i][j1].re 	+  caNoPad[i][j1].im	 * caNoPad[i][j1].im);
			caNoPadMag[i][j1]  	= sqrt(caNoPad[i1][j].re 	* caNoPad[i1][j].re 	+  caNoPad[i1][j].im  * caNoPad[i1][j].im);
			j1++;
		}
	}
}



/*************************************************************
   Over sample correlation
**************************************************************/
static void overSampleC(TrackParams *trackPar)
{
	extern fftw_complex **psFast,**psFastOver;
	extern fftw_complex **cFast,**cFastOver;
	int i,j;
	int i1,i2,j1,j2;
	fftwnd_one(trackPar->cForwardFast,cFast[0],psFast[0]);
	/*
	  Zero pad fft for over sampling
	*/
	for(i=0; i < NFAST/2; i++ ) {
		i1=NFAST*NOVER - NFAST/2 +i;
		i2=NFAST/2 + i;
		j1=NFAST*NOVER - NFAST/2;
		j2=NFAST/2;
		for(j=0; j < NFAST/2; j++ ) {
			psFastOver[i][j].re = psFast[i][j].re; 
			psFastOver[i][j].im = psFast[i][j].im;
			psFastOver[i1][j].re = psFast[i2][j].re;
			psFastOver[i1][j].im = psFast[i2][j].im;
			psFastOver[i1][j1].re = psFast[i2][j2].re;
			psFastOver[i1][j1].im = psFast[i2][j2].im;
			psFastOver[i][j1].re = psFast[i][j2].re;
			psFastOver[i][j1].im = psFast[i][j2].im;
			j1++; j2++;
		} /* End for j */
	} /* End for i */
	fftwnd_one(trackPar->cReverseFast,psFastOver[0],cFastOver[0]);
}





/******************************************************************************
    Do plans for fft's
    ************************************************************** */

static void fftCorrPlans(TrackParams *trackPar)
{
	extern fftwnd_plan aForwardIn;
	extern fftwnd_plan aReverseNoPad;
	fprintf(stderr,"1 %i %i\n", NFAST,NFAST);
	trackPar->cForwardFast = fftw2d_create_plan(NFAST, NFAST, FFTW_FORWARD,FFTW_MEASURE );
	fprintf(stderr,"3 %i %i\n",NOVER*NFAST,NOVER*NFAST);
	trackPar->cReverseFast = fftw2d_create_plan(NOVER * NFAST, NOVER*NFAST, FFTW_BACKWARD, FFTW_MEASURE);
	fprintf(stderr,"4 %i \n",trackPar->wA);
	trackPar->onedForward1 = fftw_create_plan_specific(trackPar->wA,
							   FFTW_FORWARD, FFTW_MEASURE,patch1[0],trackPar->wR,f1[0],trackPar->wR);
	trackPar->onedForward1R = fftw_create_plan_specific(trackPar->wR, FFTW_FORWARD, FFTW_MEASURE,patch1[0],1,f1[0],1);
	fprintf(stderr,"5 %i\n",trackPar->wA );
	trackPar->onedForward2 = fftw_create_plan_specific(trackPar->wA, FFTW_FORWARD, FFTW_MEASURE,patch2[0],trackPar->wR,f2[0],trackPar->wR);
	trackPar->onedForward2R = fftw_create_plan_specific(trackPar->wR, FFTW_FORWARD, FFTW_MEASURE,patch2[0],1,f2[0],1);
	fprintf(stderr,"6- %i %i\n", trackPar->wAa*OS, trackPar->wRa*OS);
	aForward = fftw2d_create_plan(trackPar->wAa*OS, trackPar->wRa*OS,  FFTW_FORWARD,FFTW_MEASURE ); /* ^^^ */
	fprintf(stderr,"8 - %i %i \n",trackPar->wAa*OS, trackPar->wRa*OS);
	aReverseNoPad = fftw2d_create_plan(trackPar->wAa*OS, trackPar->wRa*OS, FFTW_BACKWARD, FFTW_MEASURE); /* ^^^ */

	fprintf(stderr,"10 - %i %i\n", trackPar->wAa,trackPar->wRa );
	aForwardIn = fftw2d_create_plan(trackPar->wAa, trackPar->wRa,FFTW_FORWARD,FFTW_MEASURE ); /* ^^^ */
}


/*************************************************************
   Malloc space for all the global arrays/matrices and elements of trackPar
*************************************************************/

static void mallocSpace(TrackParams *trackPar)
{
	extern StrackBuf imageBuf1;
	extern StrackBuf imageBuf2;     
	extern strackComplex *cBuf1,*cBuf2;      /* Input buffers */
	extern fftw_complex **fftF1,**fftF2;
	extern fftw_complex **fftF1os,**fftF2os;
	extern fftw_complex **fftFa1,**fftFa2;
	extern fftw_complex  **psNoPad;
	extern fftw_complex **patch1, **patch2;
	extern fftw_complex **patch1in, **patch2in;
	extern fftw_complex **cNoPad;
	extern float **cFastOverMag, **cNoPadMag; 
	extern fftw_complex **cFast,**cFastOver;
	extern fftw_complex **f1, **f2;
	extern fftw_complex **img1, **img2;/* Detected images for amplitude match*/
	extern fftw_complex **img1in, **img2in;/* ^^^ Detected images for amplitude match*/
	extern fftw_complex **fftFa1os, **fftFa1Los,**fftFa2os, **fftFa2Los; /*  ^^^ */
	extern fftw_complex **psAmpNoPad;
	extern fftw_complex **caNoPad;
	extern fftw_complex **img1L, **img2L;/*images for amplitude match*/
	extern fftw_complex **psAmpNoPadL;
	extern fftw_complex **fftFa1L,**fftFa2L;      /* FFt's */
	extern fftw_complex **caNoPadL;
	extern float **caNoPadMagL;
	extern float  **caNoPadMag;     
	extern double **meanS;  
	extern double **sigmaS;
	extern double **corrResult;
	extern float **dataS,**dataR;
	fftw_complex **fftwTmp;
	strackComplex **ers1Tmp;
	int wA2,wR2;
	int i, j;
	/*
	  Input buffer.
	*/    
	imageBuf1.nr = trackPar->imageP1.nSlpR;
	imageBuf1.na = NBUFFERLINES;
	imageBuf1.firstRow = LARGEINT;
	imageBuf1.lastRow = -1;
	imageBuf2.nr = trackPar->imageP2.nSlpR;
	imageBuf2.na = NBUFFERLINES;
	imageBuf2.firstRow = LARGEINT;
	imageBuf2.lastRow = -1;
	/*
	  imageBuf1.buf=mallocfftw_complexMat(imageBuf1.na,imageBuf1.nr);
	  imageBuf2.buf=mallocfftw_complexMat(imageBuf2.na,imageBuf2.nr);
	  for(i=0; i < NBUFFERLINES; i++) 
	  for(j=0; j < imageBuf1.nr; j++) { imageBuf1.buf[i][j].re=1; imageBuf1.buf[i][j].re=0;}
	  for(i=0; i < NBUFFERLINES; i++) 
	  for(j=0; j < imageBuf2.nr; j++) { imageBuf2.buf[i][j].re=1; imageBuf2.buf[i][j].re=0;}
	*/
	if(trackPar->floatFlag ==TRUE) {
		fprintf(stderr,"Floating point input \n");
		imageBuf1.buf=(void **)mallocfftw_complexMat(imageBuf1.na,imageBuf1.nr);
		fftwTmp=(fftw_complex **)imageBuf1.buf;
		for(i=0; i < NBUFFERLINES; i++) 
			for(j=0; j < imageBuf1.nr; j++) {fftwTmp[i][j].re=1; fftwTmp[i][j].im=0;}

		imageBuf2.buf=(void**)mallocfftw_complexMat(imageBuf2.na,imageBuf2.nr);
		fftwTmp=(fftw_complex **)imageBuf2.buf;
		for(i=0; i < NBUFFERLINES; i++) 
			for(j=0; j < imageBuf1.nr; j++) {fftwTmp[i][j].re=1; fftwTmp[i][j].im=0;}

	}	else {
		fprintf(stderr,"Short int input \n");
		imageBuf1.buf=(void *)mallocErs1ComplexMat(imageBuf1.na,imageBuf1.nr);
		ers1Tmp=(strackComplex **)imageBuf1.buf;
		for(i=0; i < NBUFFERLINES; i++) 
			for(j=0; j < imageBuf1.nr; j++) {ers1Tmp[i][j].r=1; ers1Tmp[i][j].i=0;}

		imageBuf2.buf=(void*)mallocErs1ComplexMat(imageBuf2.na,imageBuf2.nr);
		ers1Tmp=(strackComplex **)imageBuf2.buf;
		for(i=0; i < NBUFFERLINES; i++) 
			for(j=0; j < imageBuf1.nr; j++) {ers1Tmp[i][j].r=1; ers1Tmp[i][j].i=0;}
	}

	wA2=(trackPar->wAa-2*trackPar->edgePadA)*OS;
	wA2=(trackPar->wRa-2*trackPar->edgePadR)*OS;

	/*
	  Patches
	*/
	patch1= mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	patch2= mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	patch1in= mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	patch2in= mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	img1= mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); /* ^^^ */
	img2= mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); /* ^^^ */
	img1in= mallocfftw_complexMat(trackPar->wAa,trackPar->wRa); /* ^^^ */
	img2in= mallocfftw_complexMat(trackPar->wAa,trackPar->wRa); /* ^^^ */
	for(i=0; i < trackPar->wAa*OS; i++) 
		for(j=0; j < trackPar->wRa*OS; j++){img1[i][j].im = 0.; img2[i][j].im = 0.;}
	/* 
	   One d strided for carrier estimation
	*/
	f1= mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	f2= mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	/*
	  FFTs
	*/
	fftF1 = mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	fftF2 = mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	fftF1os = mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	fftF2os = mallocfftw_complexMat(trackPar->wA,trackPar->wR);
	fftFa1 = mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); /* ^^^ */
	fftFa2 = mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); /* ^^^ */
	fftFa1L = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); /* ^^^ */
	fftFa2L = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); /* ^^^ */
	/* ^^^ */
	fftFa1os =  mallocfftw_complexMat(trackPar->wAa,trackPar->wRa);
	fftFa1Los = mallocfftw_complexMat(trackPar->wAa*LA,trackPar->wRa*LA);
	fftFa2os =  mallocfftw_complexMat(trackPar->wAa,trackPar->wRa);
	fftFa2Los = mallocfftw_complexMat(trackPar->wAa*LA,trackPar->wRa*LA);
	/* ^^^ */

	psNoPad = mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA);
	psFast = mallocfftw_complexMat(NFAST,NFAST);
	psFastOver = mallocfftw_complexMat(NFAST*NOVER,NFAST*NOVER);
	psAmpNoPad = mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS);
	psAmpNoPadL = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS);

	for(i=0; i < NFAST*NOVER; i++) 
		for(j=0; j < NFAST*NOVER; j++) {
			psFastOver[i][j].re=0.0; psFastOver[i][j].im=0.0;
		}
	/*    c = mallocFloatMat(trackPar->wA*NOVER,trackPar->wR*NOVER); */
	cFast = mallocfftw_complexMat(NFAST,NFAST);
	cFastOver = mallocfftw_complexMat(NFAST*NOVER,NFAST*NOVER);
	cNoPad = mallocfftw_complexMat(trackPar->wA*OSA,trackPar->wR*OSA); 
	caNoPad = mallocfftw_complexMat(trackPar->wAa*OS,trackPar->wRa*OS); 
	caNoPadL = mallocfftw_complexMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); 
	cNoPadMag = mallocFloatMat(trackPar->wA*OSA,trackPar->wR*OSA); 
	caNoPadMag = mallocFloatMat(trackPar->wAa*OS,trackPar->wRa*OS); 
	caNoPadMagL = mallocFloatMat(trackPar->wAa*LA*OS,trackPar->wRa*LA*OS); 
	cFastOverMag = mallocFloatMat(NFAST*NOVER,NFAST*NOVER);

	for(i=0; i < NFAST*NOVER; i++) 
		for(j=0; j < NFAST*NOVER; j++) {
			cFastOver[i][j].re=0.0; cFastOver[i][j].im=0.0;
			cFastOverMag[i][j]=0.0;
		}
	/* Input buffers */
	cBuf1=(strackComplex *) malloc(sizeof(strackComplex)*max(trackPar->wR,trackPar->wRa*LA));
	cBuf2=(strackComplex *) malloc(sizeof(strackComplex)*max(trackPar->wR,trackPar->wRa*LA));
	cfBuf1=(fftw_complex *) malloc(sizeof(fftw_complex)*max(trackPar->wR,trackPar->wRa*LA));
	cfBuf2=(fftw_complex *) malloc(sizeof(fftw_complex)*max(trackPar->wR,trackPar->wRa*LA));
	/* Malloc outputs.    */
	trackPar->offR =mallocFloatMat(trackPar->nA,trackPar->nR);
	trackPar->offA =mallocFloatMat(trackPar->nA,trackPar->nR);
	trackPar->corr =mallocFloatMat(trackPar->nA,trackPar->nR);
	trackPar->type =mallocByteMat(trackPar->nA,trackPar->nR);

}


/*************************************************************
   Output .dat file for matches 
*************************************************************/
static void writeDatFile(TrackParams *trackPar)
{
	FILE *fp;
	fp = fopen(trackPar->outFileD,"w");
	fprintf(fp,"%i %i %i %i %i %i\n",trackPar->rStart*trackPar->scaleFactor,trackPar->aStart*trackPar->scaleFactor,
		trackPar->nR,trackPar->nA,trackPar->deltaR*trackPar->scaleFactor,trackPar->deltaA*trackPar->scaleFactor);

}


/*************************************************************
   Read mask
*************************************************************/
static int maskValue(TrackParams *trackPar,int r1,int a1)
{
	int ia1,ir1;
	int dr,da;
	/* not mask, always try match */
	if(trackPar->maskFlag==FALSE) return(1);
	/* otherwise get mask value */
	/* fixed coords 10/29/2018 */
	dr=(trackPar->wRa/2 - trackPar->maskDat.r0)*trackPar->scaleFactor;
	da=(trackPar->wAa/2 - trackPar->maskDat.a0)*trackPar->scaleFactor;
        ia1 = (a1+da)/trackPar->maskDat.nal;
	ir1 = (r1+dr)/trackPar->maskDat.nrl;
	/* 	ia1 = (a1+trackPar->wA/2)/trackPar->maskDat.nal;	ir1 = (r1+trackPar->wR/2)/trackPar->maskDat.nrl; */
	if(ia1 < 0 || ia1 >= trackPar->maskDat.na || ir1 < 0 || ir1 >= trackPar->maskDat.nr) { 
		/*      fprintf(stderr,"ia1 ir1 %i %i\n",ia1,ir1);*/
		return(0);

	}
	return((int)(trackPar->maskDat.mask[ia1][ir1]));
}


/*************************************************************
   Routines to Malloc matrices 
*************************************************************/
static float **mallocFloatMat(int nA,int nR)
{
	float *tmp, **tmp1;
	int i;
	tmp = malloc(nR*nA*sizeof(float));
	tmp1 = (float **)malloc(nA*sizeof(float *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
}

static double **mallocDoubleMat(int nA,int nR)
{
	double *tmp, **tmp1;
	int i;
	tmp = malloc(nR*nA*sizeof(double));
	tmp1 = (double **)malloc(nA*sizeof(double *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
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

                         
static fftw_complex **mallocfftw_complexMat(int nA,int nR)
{
	fftw_complex *tmp, **tmp1;
	int i;
	tmp = malloc(nA*nR*sizeof(fftw_complex));
	tmp1 = (fftw_complex **)malloc(nA*sizeof(fftw_complex *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
}


static  strackComplex **mallocErs1ComplexMat(int nA,int nR)
{
	strackComplex *tmp, **tmp1;
	int i;
	tmp = malloc(nA*nR*sizeof(strackComplex));
	tmp1 = (strackComplex **)malloc(nA*sizeof(strackComplex *));
	for(i=0; i < nA; i++) tmp1[i]=&(tmp[i*nR]);
	return tmp1;
}





