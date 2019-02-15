#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "mosaicSource/common/common.h"
#include "math.h"
#include "stdlib.h"
#include "cRecipes/nrutil.h"
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
/* #include "ers1/landsatMosaic/landSatMosaic.h"*/
#include "cullls.h"


typedef struct intDataType {
        float x;
        float y;
} intData;

static void planeCoeffs(void *x,int i,double *afunc, int ma);

void cullLSStats(CullLSParams *cullPar)
{
	int i,j;
	int i1,j1,i2,j2;
	int ii,jj;
	unsigned long ngood;
	int good;
	double varY, varX;
	double meanY,meanX;
	double *listX,*listY,*dTmp;
	intData *data;
	unsigned long midIndex;
	float oX,oY;
	double aX[5],aY[5];
	int ma;
	int32 nx,ny;
	double chisq;
	double **u, **v, *w;    /* Tmp arrays use by svd */

	cullPar->bR=9;
	cullPar->bA=9;
	nx=cullPar->matches.nx;
	ny=cullPar->matches.ny;
	/*
	  malloc buffers
	*/
	listX=(double *) malloc(sizeof(double)*(1+(cullPar->bR+1)*(cullPar->bA +1)));
	listY=(double *) malloc(sizeof(double)*(1+(cullPar->bR+1)*(cullPar->bA +1)));
	dTmp=(double *) malloc(sizeof(double)*(1+(cullPar->bR+1)*(cullPar->bA +1)));
	for(i=0; i < (1+(cullPar->bR+1)*(cullPar->bA +1)); i++) dTmp[i]=1.0;
	data=(intData *) malloc(sizeof(intData)*(1+(cullPar->bR+1)*(cullPar->bA +1)));
	/*
	  Setup data for plane fit
	*/
	ma=3;
	u = dmatrix(1, (1+(cullPar->bR+1)*(cullPar->bA +1))+1,1,ma);
	v = dmatrix(1,ma,1,ma);
	w = dvector(1,ma);
	/*
	  Main loop over points
	*/

	for(i=0; i < ny; i++) {
		i1=max(0,i-cullPar->bA/2);
		i2=min(ny-1,i+cullPar->bA/2);

		for(j=0; j < nx; j++) {
			j1=max(0,j-cullPar->bR/2);
			j2=min(nx-1,j+cullPar->bR/2);  
			ngood=0;
			/* Save points to estimate plain, and compute mean, and list of good points */
			meanY=0.0; meanX=0.0;
			for(ii=i1; ii <= i2; ii++ ) {
				for(jj=j1; jj <= j2; jj++ ) {
					if(cullPar->matches.X[ii][jj] > (1-LARGEINT) &&
					  cullPar->matches.Y[ii][jj] > (1-LARGEINT) ) {
						ngood++;
						listX[ngood] = cullPar->matches.X[ii][jj];
						listY[ngood] =cullPar->matches.Y[ii][jj];
						meanY +=cullPar->matches.Y[ii][jj];  
						meanX +=cullPar->matches.X[ii][jj];
						data[ngood].x=(float)ii;
						data[ngood].y=(float)jj;
					}
				}
			}
			meanY /= (double)ngood; meanX /= (double)ngood;
			/* 
			   Estimate plain 
			*/
			if(ngood > 8) {
				data[0].x=meanX;
				svdfit((void *)data,listX,dTmp,ngood,aX,ma,u,v,w, 
				       &chisq, &planeCoeffs);
				data[0].x=meanY;
				svdfit((void *)data,listY,dTmp,ngood,aY,ma,u,v,w, 
				       &chisq, &planeCoeffs);
				aX[1] *=meanX; 
				aY[1] *=meanY;  
			}
			/* Compute stats after removing plain */
			ngood=0;
			meanY=0.0; meanX=0.0; varX=0.0; varY=0.0;
			for(ii=i1; ii <= i2; ii++ ) {
				for(jj=j1; jj <= j2; jj++ ) {
					if(cullPar->matches.X[ii][jj] > (1-LARGEINT) &&
					  cullPar->matches.Y[ii][jj] > (1-LARGEINT) ) {
						ngood++; 
						oX=cullPar->matches.X[ii][jj];
						oY=cullPar->matches.Y[ii][jj];
						oX-=(aX[1] + aX[2]*(float)ii + aX[3]*(float)jj);
						oY-=(aY[1] + aY[2]*(float)ii + aY[3]*(float)jj);
						meanY += oY;   meanX += oX;
						varY+=oY*oY;   varX+=oX*oX;
					} /* Endif cullPar.. */
				} /* End jj */
			} /* End ii */
			if(ngood > 8 && cullPar->matches.X[i][j] > (1-LARGEINT)) {
				cullPar->matches.sigmaY[i][j] = varY/(float)(ngood-3);
				cullPar->matches.sigmaX[i][j] = varX/(float)(ngood-3);
				cullPar->matches.sigmaY[i][j] -= meanY*meanY /
					(float)((ngood-3)*ngood);
				cullPar->matches.sigmaX[i][j] -= meanX*meanX/
                                        (float)((ngood-3)*ngood);
				cullPar->matches.sigmaY[i][j]=max(0.00001,cullPar->matches.sigmaY[i][j]);
				cullPar->matches.sigmaX[i][j]=max(0.00001,cullPar->matches.sigmaX[i][j]);
				/* Before taking square root, add variance due to quanitization error delta^2/12) */
					cullPar->matches.sigmaX[i][j] +=((1.0/NOVER)* (1.0/NOVER))/12.;
					cullPar->matches.sigmaY[i][j] += ((1.0/NOVER)* (1.0/NOVER))/12.;
				cullPar->matches.sigmaY[i][j] = sqrt((double)cullPar->matches.sigmaY[i][j]);
				cullPar->matches.sigmaX[i][j] = sqrt((double)cullPar->matches.sigmaX[i][j]);
			} else {
				cullPar->matches.sigmaY[i][j] = (float)-LARGEINT;
				cullPar->matches.sigmaX[i][j] = (float)-LARGEINT;
			}
		} /* End for j */
	} /* End for i */
}


static  void planeCoeffs(void *x,int i,double *afunc, int ma) 
{
	/*
	  Plane equation
	*/
	intData *xy; 
	xy = (intData *)x;
	/*   fprintf(stderr,"%i %f\n",i,xy[0].x);*/

	afunc[1] = xy[0].x;
	afunc[2] = (double)(xy[i].x);
	afunc[3] = (double)(xy[i].y); 
	return;
}
