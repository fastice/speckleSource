#define MXST 	75    /* maximum number of orbital state vectors   
                        in the processing parameter file */


typedef struct{double x,y,z;} VEC; /* vector structure cartesian (X,Y,Z)*/

typedef struct{			  /* state vector structure */
   VEC pos;			  /* position vector */
   VEC vel;			  /* velocity vector */
} STATE;

typedef struct {
   int hour;
   int minute;
   double second;
} timehms;
                   
typedef struct {
  char title[512];	/* ascii text description of scene and processing 
                           parameters */
  char date[16];	/* date data acquired/generated in decimal numbers 
                           (DD MM YYYY) */
  char asciitime[24];	/* time of first record in SAR signal data file 
                           (HH MM SS.SSSS) */
  timehms time;         /* Time in hours minutes seconds (int,int,fp) */
  char pol_ch[16];	/* polarization/channel to be processed out of the set 
			   {HH, HV, VH, VV, CH1, CH2} */
  double lat,lon;	/* latitude and longitude of scene center 
                            (decimal degrees) */
  double track;		/* track angle of radar platform (decimal degrees) */
  double alt;		/* average altitude above geoid of the radar 
                           platform (m) */
  double terra_alt;	/* average reference height of terrain in the scene
                           above the geoid (m) */
  VEC pos;		/* position of platform at center of presum data in
                           (X,Y,Z/T,C,N) (m) */
  VEC vel;		/* velocity vector (X,Y,Z/T,C,N) along track (m/s) 
                           at center of presum data */
  VEC acc;		/* acceleration vector (X,Y,Z/T,C,N) along track 
                           (m/s^2) at center of presum data */
  double prf;		/* radar pulse repetition frequency PRF (Hz) */
  double fdp[4];	/* doppler centroid polynomial coefficients as a
                           function of range
                           fd=fdp[0]+fdp[1]*r+fdp[2]*r**2+fdp[3]*r**3 (Hz)/m*/
  double td;		/* time delay between transmission of pulse and 
                           first sample of echo (seconds) */ 
  int nr_offset;	/* offset in (samples) to first range sample to process,
			   each I/Q pair counts as 1 sample */	
  double r0,r1,r2;	/* raw SAR data near, center, far slant ranges (m) */
  double ir0,ir1,ir2;	/* output image near, center and far range (m) */
  double rpixsp,ran_res; /* slant-range pixel spacing, slant-range resolution 
                           (m) */
  double prfrac;	/* fraction of doppler bandwidth to process 
                           (0.0 -- 1.0) */ 
  int nprs_az;		/* initial azimuth presumming factor 
                           (4 for a 4:1 presum, for example) */
  int loff;		/* offset in echoes from the start of SAR signal 
                           data file to the first presum output record */
  int nl;		/* # echoes to process/simulate */
  int nrfft;		/* # of range samples to process (set equal to the 
                           range FFT size, power of 2) */
  int nsub;		/* # sub-apertures to divide the presummed doppler
                           spectrum */
  int nlr,nlaz;		/* # range looks, # azimuth looks/sub-aperture */
  double azoff;		/* along track azimuth offset of the first image line 
                           from the start of the presummed data (m) */
  double azimsp, azimres; /* azimuth image pixel spacing, azimuth image
                             resolution (m) */
  int nrs, nazs;          /* image width in range pix, length in azimuth pix */
  int nstate;		  /* number of state vectors */
  double t_state;	  /* UTC time (sec) since start of day for 
                             first state vector */
  double tis;		  /* time interval between state vectors (s) */
  double sTimes[MXST+1];
  STATE state[MXST+1];	  /* maximum of MXST state vectors (X,Y,Z) 
                             convential terrestrial system (CTS) */
  double lambda;          /* Wavelength */  
  double lookDir;         /* +1.0 for right looking, -1 for left looking */
} PROC_PAR;

   void parsePar(char *parFile, PROC_PAR *procp);

  void glatlon( PROC_PAR *par, double ***lat1,double ***lon1,int ma, int mr,
 double deltaT);

    void smlocate(double rxs,  double rys,  double rzs,
                  double rvxs, double rvys, double rvzs,
                  double rsl, double fd,double *lat,double *lon, PROC_PAR *par);

    void centerLL( PROC_PAR *par, double *lat, double *lon, double deltaT);
    void correctTime(PROC_PAR *par, double *squint, double noffset,
                     double *tskew,double *toffset, int squintTime);
