#include "stdio.h"
#include"string.h"
#include "clib/standard.h"
#include "cullst.h"
#include "math.h"
#include "stdlib.h"
#include "mosaicSource/common/common.h"
#include "cRecipes/nrutil.h"

typedef struct intDataType {
        float x;
        float y;
} intData;

 static void planeCoeffs(void *x,int i,double *afunc, int ma);

    void cullStats(CullParams *cullPar)
{
    int i,j;
    int i1,j1,i2,j2;
    int ii,jj;
    unsigned long ngood;
    int good;
    double varA, varR;
    double sigmaA,sigmaR;
    double meanA,meanR;
    double *listR,*listA,*dTmp;
    intData *data;
    unsigned long midIndex;
    float oR,oA;
    double aR[5],aA[5];
    int ma;
    double chisq;
    double **u, **v, *w;    /* Tmp arrays use by svd */

    cullPar->bR=9;
    cullPar->bA=9;
    listR=(double *)
        malloc(sizeof(double)*(1+(cullPar->bR+1)*(cullPar->bA +1)));
    listA=(double *)
        malloc(sizeof(double)*(1+(cullPar->bR+1)*(cullPar->bA +1)));
    dTmp=(double *)
        malloc(sizeof(double)*(1+(cullPar->bR+1)*(cullPar->bA +1)));
    for(i=0; i < (1+(cullPar->bR+1)*(cullPar->bA +1)); i++) dTmp[i]=1.0;
    data=(intData *)
        malloc(sizeof(intData)*(1+(cullPar->bR+1)*(cullPar->bA +1)));

     ma=3;
     u = dmatrix(1, (1+(cullPar->bR+1)*(cullPar->bA +1))+1,1,ma);
     v = dmatrix(1,ma,1,ma);
     w = dvector(1,ma);


    for(i=0; i < cullPar->nA; i++) {
        i1=max(0,i-cullPar->bA/2);
        i2=min(cullPar->nA-1,i+cullPar->bA/2);
        for(j=0; j < cullPar->nR; j++) {
            j1=max(0,j-cullPar->bR/2);
            j2=min(cullPar->nR-1,j+cullPar->bR/2);  
            ngood=0;
            /* Save points to estimate plain */
            meanA=0.0; meanR=0.0;
            for(ii=i1; ii <= i2; ii++ ) {
                for(jj=j1; jj <= j2; jj++ ) {
                    if(cullPar->offR[ii][jj] > (1-LARGEINT) &&
                       cullPar->offA[ii][jj] > (1-LARGEINT) ) {
                        ngood++;
                        listR[ngood] = cullPar->offR[ii][jj];
                        listA[ngood] = cullPar->offA[ii][jj];
                        meanA += cullPar->offA[ii][jj];  
                        meanR +=cullPar->offR[ii][jj];
                        data[ngood].x=(float)ii;
                        data[ngood].y=(float)jj;
                    }
                }
            }
            meanA /= (double)ngood; meanR /= (double)ngood;
            /* Estimate plain */
            if(ngood > 8) {
	      data[0].x=meanR;
                svdfit((void *)data,listR,dTmp,ngood,aR,ma,u,v,w, 
                   &chisq, &planeCoeffs);
	      data[0].x=meanA;
                svdfit((void *)data,listA,dTmp,ngood,aA,ma,u,v,w, 
		&chisq, &planeCoeffs);
	        aR[1] *=meanR; 
	        aA[1] *=meanA;  
            }
	    /* Compute stats after removing plain */
            ngood=0;
            meanA=0.0; meanR=0.0; varR=0.0; varA=0.0;
            for(ii=i1; ii <= i2; ii++ ) {
                for(jj=j1; jj <= j2; jj++ ) {
                    if(cullPar->offR[ii][jj] > (1-LARGEINT) &&
                       cullPar->offA[ii][jj] > (1-LARGEINT) ) {
                        ngood++; 
                        oR=cullPar->offR[ii][jj];
                        oA=cullPar->offA[ii][jj];
                        oR-=(aR[1] + aR[2]*(float)ii + aR[3]*(float)jj);
                        oA-=(aA[1] + aA[2]*(float)ii + aA[3]*(float)jj);
                        meanA += oA;   meanR += oR;
                        varA+=oA*oA;   varR+=oR*oR;
                    } /* Endif cullPar.. */
	          } /* End jj */
            } /* End ii */
 	    if(ngood > 8 && cullPar->offR[i][j] > (1-LARGEINT)) {
                cullPar->sigmaA[i][j] = varA/(float)(ngood-3);
                cullPar->sigmaR[i][j] = varR/(float)(ngood-3);
                cullPar->sigmaA[i][j] -= meanA*meanA /
                                         (float)((ngood-3)*ngood);
                cullPar->sigmaR[i][j] -= meanR*meanR/
                                        (float)((ngood-3)*ngood);
		cullPar->sigmaA[i][j]=max(0.00001,cullPar->sigmaA[i][j]);
		cullPar->sigmaR[i][j]=max(0.00001,cullPar->sigmaR[i][j]);
                cullPar->sigmaA[i][j] = sqrt((double)cullPar->sigmaA[i][j]);
                cullPar->sigmaR[i][j] = sqrt((double)cullPar->sigmaR[i][j]);
	    } else {
                cullPar->sigmaA[i][j] = (float)-LARGEINT;
                cullPar->sigmaR[i][j] = (float)-LARGEINT;
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
