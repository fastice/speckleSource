#include "stdio.h"
#include"string.h"
#include "mosaicSource/common/common.h"
#include "clib/standard.h"
#include "math.h"
#include "geotiff/xtiffio.h"  /* for TIFF */
#include "geotiff/geotiffio.h" /* for GeoTIFF */
#include "landsatSource64/Lstrack/lstrack.h"
#include "landsatSource64/Lsfit/lsfit.h"
/* #include "ers1/landsatMosaic/landSatMosaic.h"*/
#include "cullls.h"

    void cullLSSmooth(CullLSParams *cullPar)
{
    int i,j;
    int i1,j1,i2,j2;
    int ii,jj;
    float nlooks[4];
    unsigned long ngood;
    double meanY,meanX;
    float nLooksEff;
 int32 nx,ny;
    nlooks[0]=1.0;
    nlooks[1]=1.0;
    nlooks[2]=1.25;
    nlooks[3]=1.015;
    nx=cullPar->matches.nx;
    ny=cullPar->matches.ny;
    fprintf(stderr,"\n**** smoothing data sr=%i  sa=%i ****\n\n",  cullPar->sX,cullPar->sY);

    for(i=0; i < ny; i++) {
        i1=max(0,i-cullPar->sY/2);
        i2=min(ny-1,i+cullPar->sY/2);
/*        if(cullPar->sY == 1) {i1=1; i2=1;}*/
        for(j=0; j < nx; j++) {
            j1=max(0,j-cullPar->sX/2);
            j2=min(nx-1,j+cullPar->sX/2);
/*            if(cullPar->sX == 1) {j1=1; j2=1;}*/
            ngood=0;
            nLooksEff=0.0;
            meanY=0.0; meanX=0.0; 
            for(ii=i1; ii <= i2; ii++ ) {
                for(jj=j1; jj <= j2; jj++ ) {
                    if(cullPar->matches.X[ii][jj] > (1-LARGEINT) &&  cullPar->matches.Y[ii][jj] > (1-LARGEINT) ) {
                        ngood++; 
			nLooksEff++; /* NEED TO UPDATE THIS */
/*                        nLooksEff += nlooks[(int)(cullPar->type[ii][jj])];*/
                        meanY += cullPar->matches.Y[ii][jj];
                        meanX += cullPar->matches.X[ii][jj];

                    } /* Endif cullPar.. */
	          } /* End jj */
            } /* End ii */
/* added exception for sX=1 or sY=1 */
 	    if(((ngood > 1)||(cullPar->sY==1)||(cullPar->sX==1)) && cullPar->matches.X[i][j] > (1-LARGEINT) &&  cullPar->matches.sigmaY[i][j] > (1-LARGEINT)) {
   	        nLooksEff=max(nLooksEff,1.);
                cullPar->matches.sigmaY[i][j] /= sqrt((double)nLooksEff);
                cullPar->matches.sigmaX[i][j] /= sqrt((double)nLooksEff);
                cullPar->XS[i][j] = meanX/(float)ngood;
                cullPar->YS[i][j] = meanY/(float)ngood;
	    } else {
                cullPar->XS[i][j] = (float)-LARGEINT;
                cullPar->YS[i][j] = (float)-LARGEINT;
            }
        } /* End for j */
    } /* End for i */
}


