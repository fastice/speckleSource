#include "speckleSource/Strack/strack.h"
#include "cRecipes/nrutil.h"
/*
  sl pair structure with int values
*/
typedef struct slIntPairType {
   int s;     /* s Value */
   int l;     /* l Value */
} slIntPair;

/*
   Parameters for match
*/

typedef struct correlationResultType {
   slIntPair loc1;
   slIntPair loc2;
   slIntPair delta;
   double meanS;
   double meanR;
   double varS;
   double varR;
   double corrRS;
   double corrPt;
   double cov[4];
   double snr;
 } correlationResult;
