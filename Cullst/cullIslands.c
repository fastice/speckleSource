#include "stdio.h"
#include "string.h"
#include "clib/standard.h"
#include "gdalIO/gdalIO/grimpgdal.h"
#include "cullst.h"
#include "utilsSource/intFloat/intFloat.h"
#include <stdlib.h>
#include "math.h"
#include "unwrapSource/unWrap/unwrapInclude.h"

static void screenIslands(hole *holes, int32_t maxLabel, int32_t islandThresh);

void cullIslands(CullParams *cullPar)
{
    unwrapImageStructure uwImage; /* This is to fake out for using
                                unwrap image connect components */
    hole *holes;
    int32_t ncull, npixcull;
    int32_t x, y;   /* Current location in hole */
    int32_t xp, yp; /* Hole point where last update occured */
    double a[8];
    char **mask;         /* Mask used for mark pixels for cc alg */
    double *zTmp, *dTmp; /* Data points use in interpolation for curr pt */
    double **u, **v, *w; /* Tmp arrays use by svd */
    double dh;           /* Distance between consecutive points in hole */
    int32_t **labels;    /* connect components labeled output */
    int32_t nr, na;      /* Size in range and azimuth */
    int32_t ma;          /* Number of parameters interpolater estimate */
    int32_t i, j, k;     /* Loop variables */
    int32_t width, height;
    int32_t maxLabel;        /* Largest labels from cc */
    int32_t *buf;            /* Temp variable used for mem allocation */
    int32_t lList[25], nLab; /* List of labels for border point */
    int32_t nPts;
    float minD;
    int32_t skipValue;
    int32_t maxBorders;   /* Maximum border points over all holes */
    int32_t skip, update; /* Flags used in interpolation loop */
    FILE *fp;
    /*
       Init  image stuff
   */
    na = cullPar->nA;
    nr = cullPar->nR;
    fprintf(stderr, "Culling Islands - %i %i %i\n", cullPar->islandThresh, na, nr);
    /*
       Mask image
    */
    mask = (char **)malloc(na * sizeof(char *));
    for (i = 0; i < na; i++)
        mask[i] = (char *)malloc(nr * sizeof(char));
    for (i = 0; i < na; i++)
        for (j = 0; j < nr; j++)
        {
            mask[i][j] = VALIDPIXEL;
            if (cullPar->offRS[i][j] >= (float)(-LARGEINT + 1000))
                mask[i][j] += REGION;
        }
    /*
        Label image
    */
    labels = (int32_t **)malloc(na * sizeof(int32_t *));
    buf = malloc(na * nr * sizeof(int));
    for (i = 0; i < na; i++)
        labels[i] = (int32_t *)&(buf[i * nr]);
    /*
       Copy to uw structure for connect components.
    */
    uwImage.phase = cullPar->offRS;
    uwImage.rangeSize = nr;
    uwImage.azimuthSize = na;
    uwImage.bCuts = mask;
    uwImage.labels = labels;
    /*
       Run connected components
    */
    fprintf(stderr, "Labeling regions\n");
    labelRegions(&uwImage);
    /*
    fp = fopen("junk.labels","w");
    for(i=0; i < na; i++) fwriteBS(uwImage.labels[i],sizeof(int32_t),nr,fp,INT32FLAG); fclose(fp);
    */
    /*
        Pass 1: Get max label
    */
    fprintf(stderr, "Pass 1 - checking labels %i %i \n", nr, na);
    maxLabel = 0;
    for (i = 0; i < na; i++)
        for (j = 0; j < nr; j++)
            if (labels[i][j] < 0.99 * LARGEINT)
                maxLabel = max(maxLabel, labels[i][j]);
    fprintf(stderr, "maxLabel = %i \n", maxLabel);
    holes = (hole *)malloc(sizeof(hole) * (maxLabel + 1));
    for (i = 1; i <= maxLabel; i++)
    {
        holes[i].nH = 0;
        holes[i].nB = 0;
        holes[i].cH = 0;
        holes[i].cB = 0;
    }
    /*
       Pass 2: Count holes and border pixels
    */
    fprintf(stderr, "Pass 2 - counting holes and borders\n");
    for (i = 0; i < na; i++)
        for (j = 0; j < nr; j++)
        {
            if (labels[i][j] < LARGEINT && validPixel(mask[i][j]) == TRUE)
            {
                (holes[labels[i][j]].nH)++;
            }
        }
    /*
       Now allocate space for arrays
    */
    maxBorders = 0;
    for (i = 1; i <= maxLabel; i++)
    {
        if (holes[i].nH > 0)
        {
            holes[i].xb = (int32_t *)malloc(holes[i].nB * sizeof(int));
            holes[i].yb = (int32_t *)malloc(holes[i].nB * sizeof(int));
            holes[i].xh = (int32_t *)malloc(holes[i].nH * sizeof(int));
            holes[i].yh = (int32_t *)malloc(holes[i].nH * sizeof(int));
        }
    }
    /*
       Pass 3: Save locations of hole and border pixels
    */
    fprintf(stderr, "Pass 3 -  saving holes and borders\n");
    for (i = 0; i < na; i++)
    {
        for (j = 0; j < nr; j++)
        {
            if (labels[i][j] < LARGEINT && validPixel(mask[i][j]) == TRUE)
            {
                holes[labels[i][j]].xh[holes[labels[i][j]].cH] = j;
                holes[labels[i][j]].yh[holes[labels[i][j]].cH] = i;
                (holes[labels[i][j]].cH)++;
            }
        }
    }
    /*
        screen Holes
    */
    screenIslands(holes, maxLabel, cullPar->islandThresh);
    /*
       cull Islands
    */
    ncull = 0;
    npixcull = 0;
    for (i = 1; i < maxLabel; i++)
    {
        if (holes[i].nH > 0)
        {
            ncull++;
            for (j = 0; j < holes[i].nH; j++)
            {
                npixcull++;
                /* allow intfloat to use */
                if (cullPar->offAS != NULL)
                    cullPar->offAS[holes[i].yh[j]][holes[i].xh[j]] = (float)-LARGEINT;
                cullPar->offRS[holes[i].yh[j]][holes[i].xh[j]] = (float)-LARGEINT;
            }
        }
    }
    fprintf(stderr, "Killed %i islands with a total of %i pixels\n", ncull, npixcull);
}

/*
   Use various criteria to avoid fixing holes
*/
static void screenIslands(hole *holes, int32_t maxLabel, int32_t islandThresh)
{
    int32_t i, j;
    int32_t minX, maxX, minY, maxY;
    float ratio;

    for (i = 1; i < maxLabel; i++)
    {

        minX = LARGEINT;
        maxX = -LARGEINT;
        minY = LARGEINT;
        maxY = -LARGEINT;
        for (j = 0; j < holes[i].nH; j++)
        {
            minX = min(minX, holes[i].xh[j]);
            maxX = max(maxX, holes[i].xh[j]);
            minY = min(minY, holes[i].yh[j]);
            maxY = max(maxY, holes[i].yh[j]);
        }
        if ((maxX - minX) > islandThresh || (maxY - minY) > islandThresh)
        {
            holes[i].nH = 0;
        }
    }
}
