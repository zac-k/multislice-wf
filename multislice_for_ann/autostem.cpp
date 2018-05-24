/*      *** autostem.cpp ***

------------------------------------------------------------------------
Copyright 1998-2012 Earl J. Kirkland


This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------- NO WARRANTY ------------------
THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
IN NO EVENT SHALL THE AUTHOR BE LIABLE
FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
ANY OTHER PROGRAM). 

-----------------------------------------------------------------------------
   ANSI-C and TIFF version
  this version uses FFTW 3 (net about a factor of 2X faster)
  see:   www.fftw.org

  to output the position averaged CBED pattern include the filename
  after the program name on the command line (in 2D image mode only)

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -fopenmp -o autostem autostem.cpp slicelob.o
                      floatTIFF.o  -lfftw3f

  Calculate STEM images and line scans from a non-periodic
  distribution of (x,y,z) atomic coordinates using repetitive multislice
  
  The transmission function for each slice can take a lot of
  computer time.  This program propagates the STEM probe for many
  beam position through the specimen at the same time to avoid recalculating
  the specimen transmission function for each probe position. This
  requires a huge amount of memory.  In 1D line scan mode a 2D probe
  wave function is stored for all positions.  In 2D image mode the 2D
  probe wave functions for a whole line are stored simulataneously.

  this file is formatted for a tab size of 4 characters

  multithreaded code using openMP may be turned on/off with
  symbol USE_OPENMP (ignore undefined pragma warning if off)

  Ref:
  [1] Zhiheng Yu, Philip Batson, John Silcox, "Artifacts in Aberration
      Corrected ADF-STEM Imaging", Ultramicroscopy 96 (2003) p.275-284
  [2] J. M. LeBeau et al, "Position averaged convergent beam electron
      diff.: Theory and Applications", Ultramic. 110 (2010) p.118-125

  started from stemslic.c and autoslic.c 19-jul-1998 E. Kirkland
  convert to simultaneous transmission of many probes at the same
    time to reuse transmission functions  28-jul-1998 ejk
  finished 2D image mode (1D mode works OK) 29-jul-1998 ejk
  last modified 29-jul-1998 ejk
  add Mote-Carlo integration of source size 24-aug-1998 ejk
  fixed typo in random dither in y direction of probe position
    1-feb-1999 ejk
  fixed error in nprobes in image mode (should be nyout
      but it was nxout- at top question) 3-july-2001 ejk
  update memory allocation routines,
     change void main() to int main() for better portability,
     and add 5th order spherical aberration  3-mar-2004 ejk
  start modification to keep multiple thickness's 21-jul-2005 ejk
      - in working form 26-jul-2005 ejk
  put in faster sorting routine sortByZ() 6-feb-2006 ejk
  add some openMP multithreading 23-may-2006 ejk
  move openMP stuff into conditional compilation
     so I only need one version of the code,
     and move sortbyz() to sliclib 23-aug-2006 ejk
  add periodic() to put probe position back into the supercell
     and recal. x,y position without source size wobble
     for output in 1D mode 31-aug-2006 ejk
  change propagation range to be whole unit cell not just
     range of atoms to treat sparsely populated spec.
     better 14-sep-2006 ejk
  change range to start at -0.25*deltaz 26-sep-2006 ejk
  minor cosmetic changes 11-jul-2007 ejk
  minor fixes for openMP stuff 27-feb-2008 ejk
  move vzatomLUT() to slicelib.c and reformat for TAB size of 4 char
               11-jun-2008 ejk
  convert to GPL 3-jul-2008 ejk
  fix bug in multithreading (add more private var.) 
      5-nov-2008 ejk
  add confocal mode 31-jul-2009 E. Kirkland
  offset confocal collector by zslice (same ref as obj) 4-aug-2009 ejk
  fix param listing in confocal mode 3,6-dec-2009 ejk
  get return value of scanf() to remove warnings from gcc 4.4 
         10-feb-2010 ejk
  start conversion to FFTW  10-feb-2010 ejk
     in working form 16-feb-2010 ejk
  move vz LUT init (for openMP) to top 22-feb-2010 ejk
  fix sign convention in FFTW 21-mar-2010 ejk
  update comments 4-apr-2010 ejk
  add multipole aberrations to probe (but not collector yet)
       9-may-2011 ejk
  change detector max test to < from <= so adding many
     ADF detectors together will work 6-jul-2011 ejk
  add trap for undersampling the probe+aperture 3-mar-2012 ejk
  start conversion to floatTIFF and C++ 24-may-2012 ejk
  add option to save position averaged CBED 30-jun-2012 ejk
*/

#include <cstdio>  /* ANSI C libraries used */
#include <cstdlib>
#include <cstring>
#include <cmath> 
#include <ctime>

#include "fftw3.h"        /* FFT routines from FFTW 3*/
#include "slicelib.hpp"   /* misc. routines for multislice */
#include "floatTIFF.hpp"  /* file I/O routines in TIFF format */

#define MANY_ABERR      /*  define to include many aberrations */

//
#define USE_OPENMP      /* define to use openMP */

#ifdef USE_OPENMP
#include <omp.h>
/*  get wall time for benchmarking openMP */
#define walltim() ( omp_get_wtime() )
double walltimer;
#endif

const float BW= (2.0F/3.0F);   /* antialiasing bandwidth limit factor */
const double ABERR= 1.0e-5;    /* max error for a,b */

const int NCMAX=   512; /* max number of characters per line */

const int ADF=      0;  /* modes of collector */
const int CONFOCAL= 1;

const int TRUE=    1;
const int FALSE=   0;

const int NZMAX=   103;    /* max atomic number Z */
const int NRMAX=   100;    /* number of in look-up-table in vzatomLUT */
const double RMIN=    0.01;   /* r (in Ang) range of LUT for vzatomLUT() */
const double RMAX=    5.0;

/* global specimen and probe parameters to share */

fftwf_complex **probe;          /* complex probe wave functions */
fftwf_complex *trans;           /* complex transmission functions */
fftwf_plan  planPf, planPi;   /* FFTW plans forward/inverse for probe,  */
fftwf_plan  planTf, planTi;   /*  transmission functions */
float *propxr, *propxi;   /* complex propagator vs x */
float *propyr, *propyi;   /* complex propagator vs y */
float zmin, zmax;
float **pacbedPix;        /*  to save position averaged CBED */

int nx, ny, nxprobe, nyprobe, nslice;
int natom;
float *kx, *ky, *kx2, *ky2, *kxp, *kyp, *kxp2, *kyp2;
float  **rmin, **rmax, *xp, *yp;
float *xa, *ya, *za, *occ, *wobble;
float *xa2, *ya2, *za2, *occ2;
int *Znum, *Znum2;
double wavlen, k2maxp, Cs3,Cs5, df,apert1, apert2, pi, keV;
float ax, by, cz;                   /*  specimen dimensions */
double *almin, *almax, *k2max, *k2min, deltaz;
long nbeamt;

floatTIFF myFile;

/* define functions at end of this file (i.e. so main can be 1st) */
double periodic( double pos, double size );
void STEMsignals( double x[], double y[], int npos, float p[],
         int multiMode, double ***detect, int ndetect,
         double ThickSave[], int nThick, double sum[] );
void invert2D( float** pix, long nx, long ny );    /*   for CBED pix */

/* spline interpolation coeff. */
int splineInit=0, *nspline;
double  *splinx, **spliny, **splinb, **splinc, **splind;

void trlayer( const float x[], const float y[], const float occ[],
    const int Znum[], const int natom, 
    const float ax, const float by, const float kev,
    fftwf_complex *trans, const long nx, const long ny,
    double *phirms, long *nbeams, const float k2max );

/* extra globals for confocal mode */
int doConfocal;
int *collectorMode;
float dfa2C, dfa2phiC, dfa3C, dfa3phiC;   /* astigmatism parameters */
double *collectMin, *collectMax;
double Cs3C, Cs5C, dfC, apert1C, apert2C; /* aberrations of collector lens */

int main( int argc, char *argv[ ] )
{
    char filein[NCMAX], fileout[NCMAX], fileoutpre[NCMAX],
        description[NCMAX], cmode, cline[NCMAX], pacbedFile[NCMAX];
  
    const char version[] = "30-jun-2012 (ejk)";

    int ix, iy, i, idetect, nout, nxout, nyout,
        ncellx, ncelly, ncellz, iwobble, nwobble,
        ndetect, nprobes, ip, nThick, it, ns,
        done, status, multiMode;

    int l1d, lwobble, lxzimage, NPARAM;
    int lpacbed, ixo,iyo, ix2,iy2, nx1,nx2, ny1,ny2, nxout2,nyout2;

    long nbeamp, nbeampo;
    long  ltime;
    unsigned long iseed;

    float *param, ***pixr, temp, pmin, pmax, aimin, aimax;
    float wmin, wmax, xmin,xmax, ymin, ymax, temperature;
    float prr, pri, rmin0, rmax0, scalef;
    double dxp, dyp;

    double scale, *x, *y, sum, *sums, w, ***detect, ***detect2,
       tctx, tcty, xi,xf, yi,yf, dx, dy, totmin, totmax,
       ctiltx, ctilty, timer, sourcesize, sourceFWHM, *ThickSave,
       vz, rsq, k2maxa, k2maxb, k2;

    FILE *fp;

/* start by announcing version etc */

    printf("autostem (ADF,confocal) version dated %s\n", version );
    printf("Copyright (C) 1998-2012 Earl J. Kirkland\n" );
    printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
            " under the GNU general public license\n\n" );

    printf( "Calculate STEM images using FFTW\n");
#ifdef USE_OPENMP
    printf( "and multithreaded using openMP\n");
#endif
    printf( "\n" );

/*    get option to save position averaged CBED pattern */
    if( argc >= 2 ) {
        strncpy( pacbedFile,  argv[1], NCMAX ); 
        printf( "calculate 2D position averaged CBED pattern (arb. units) in file %s\n",
                pacbedFile );
        lpacbed = TRUE;
    } else {
        printf( "To calculate 2D pos. aver. CBED, include a file name on command line.\n");
        lpacbed = FALSE;
    }

/*  get simulation options */

    pi = 4.0 * atan( 1.0 );
    NPARAM = myFile.maxParam();
    param = (float*) malloc1D( NPARAM, sizeof(float), "param" );
    for( i=0; i<NPARAM; i++) param[i] = 0.0F;

    printf("Name of file with input atomic "
           "potential in x,y,z format:\n");
    ns = scanf("%500s", filein );

    printf("Replicate unit cell by NCELLX,NCELLY,NCELLZ :\n");
    ns = scanf("%d %d %d", &ncellx, &ncelly, &ncellz);
    if( ncellx < 1 ) ncellx = 1;
    if( ncelly < 1 ) ncelly = 1;
    if( ncellz < 1 ) ncellz = 1;

/*  get more parameter etc */

    printf("STEM probe parameters, V0(kv), Cs3(mm), Cs5(mm),"
           " df(Angstroms), apert1,2(mrad) :\n");
    ns = scanf("%lg %lg %lg %lg %lg %lg",
          &keV, &Cs3, &Cs5, &df, &apert1, &apert2);
    param[pDEFOCUS] = (float) df;
    param[pCS]  = (float) ( Cs3*1.0e7 );
    param[pCS5] = (float) ( Cs5*1.0e7 );

#ifdef MANY_ABERR
    /*   get higher order aberrations if necessary */
    printf("type higher order aber. name (as C32a, etc.) followed\n"
        " by a value in mm. (END to end)\n");
    done = multiMode = 0;
    do{
        ns = scanf( "%20s", cline );
        if( strstr( cline, "END" ) != NULL ) {
            done = 1;
        } else {
            ns = scanf( "%lg", &vz );
            status = readCnm( cline, param, vz );        
            if( status < 0 ) {
                printf( "unrecognized aberration, exit...\n");
                exit( EXIT_SUCCESS );
            } else multiMode = 1;
        }
    } while( !done );

#endif

    wavlen = wavelength( keV );
    printf("wavelength = %f Angstroms\n", wavlen);
    if( apert1 > apert2 ) {
       printf("Bad probe aperture specification.\n");
       printf("apert1 must be less than apert2.\n");
       printf("apert1=%f, apert2 = %f\n", apert1, apert2);
       exit( 0 );
    }

    /*  FFTW is not restricted to 2^n sizes so don't need test anymore */
    printf("Size of specimen transmission function"
          " Nx,Ny in pixels : \n");
    ns = scanf("%d %d", &nx, &ny);

    printf("Size of probe wave function"
          " Nx,Ny in pixels : \n");
    ns = scanf("%d %d", &nxprobe, &nyprobe);

    printf("Crystal tilt x,y in mrad. :\n");
    ns = scanf( "%lf %lf", &ctiltx, &ctilty );
    ctiltx = ctiltx * 0.001;
    ctilty = ctilty * 0.001;

    l1d = askYN("Do you want to calculate a 1D line scan");

    if( l1d == 1 ) {
        lxzimage = askYN("Do you want to save all depth information as xz image");
        nThick = 1;
    } else {
        do { printf("Number of thickness levels to save, including"
                " the end(>=1):\n");
             ns = scanf( "%d", &nThick );
        } while (nThick <= 0);
        ThickSave = (double*) malloc1D( nThick, sizeof(double), "ThickSave");
        if( nThick > 1 ) {
           printf( "type thickness (in Ang.) of %d intermediate layers"
                " :\n", (nThick-1) ); 
           for( it=0; it<(nThick-1); it++) ns = scanf( "%lf", &ThickSave[it] );
        }
    }

    printf("File name prefix to get output of STEM multislice result "
        "(no extension):\n");
    ns = scanf("%500s", fileoutpre);

    do {    printf("Number of detector geometries (>=1):\n");
            ns = scanf( "%d", &ndetect );
    } while (ndetect <= 0);

    almin  = (double*)  malloc1D( ndetect, sizeof(double), "almin" );
    almax  = (double*)  malloc1D( ndetect, sizeof(double), "almax" );
    collectorMode = (int*) malloc1D( ndetect, sizeof(int), "collectorMode" );

    doConfocal = FALSE;

    for( idetect=0; idetect<ndetect; idetect++) {
        printf("Detector %3d, type: min max angles(mrad)"
                " or radius(Ang.) \n followed by m or A\n", idetect+1);
        ns = scanf("%lg %lg %c",
                &almin[idetect], &almax[idetect], &cmode );
        if( (cmode == 'm') || (cmode=='M') ) {
                collectorMode[idetect] = ADF;
                printf( "normal ADF detector\n");
        } else if( (cmode == 'a') || (cmode=='A') ) {
                collectorMode[idetect] = CONFOCAL;
                printf( "confocal detector\n");
                doConfocal = TRUE;
        } else {
                printf( "unrecognized collector mode = %c\n", cmode);
                exit( 0 );
        }
    }  /* end for(idetect=.. */

    if( doConfocal == TRUE ) {
        printf("Collector lens parameters, Cs3(mm), Cs5(mm),"
                 " df(Angstroms), apert1,2(mrad) :\n");
        ns = scanf("%lg %lg %lg %lg %lg",
                         &Cs3C, &Cs5C, &dfC, &apert1C, &apert2C);
        printf( "Magnitude and angle of 2-fold astigmatism"
                        " (in Ang. and degrees):\n");
        ns = scanf( "%f %f", &dfa2C, &dfa2phiC);
        dfa2phiC = (float) (dfa2phiC * pi /180.0F);
        printf( "Magnitude and angle of 3-fold astigmatism"
                " (in Ang. and degrees):\n");
        ns = scanf( "%f %f", &dfa3C, &dfa3phiC);
        dfa3phiC = (float) (dfa3phiC * pi /180.0F);

        if( apert1C > apert2C ) {
                printf("Bad collector aperture specification.\n");
                printf("apert1 must be less than apert2.\n");
                printf("apert1=%f, apert2 = %f\n", apert1C, apert2C);
                exit( 0 );
        }
    }  /* end if( doConfocal...  */

    if( l1d == 1 ) {
        printf("xi, xf, yi, yf, nout :\n");
        ns = scanf("%lg %lg %lg %lg %d", &xi, &xf, &yi, &yf, &nout);
        nprobes = nout;
    }else {
        printf("xi,xf,yi,yf, nxout,nyout :\n");
        ns = scanf("%lg %lg %lg %lg %d %d",
             &xi, &xf, &yi, &yf, &nxout, &nyout);
        nprobes = nyout;
    }

    /*  remember that the slice thickness must be > atom size
        to use projected atomic potential */
    printf("Slice thickness (in Angstroms):\n");
    ns = scanf("%lf", &deltaz );
    if( deltaz < 1.0 ) {
        printf("WARNING: this slice thickness is probably too thin"
            " for autostem to work properly.\n");
    }

    lwobble = askYN("Do you want to include thermal vibrations");
    if( lwobble == 1 ) {
        printf( "Type the temperature in degrees K:\n");
        ns = scanf( "%g", &temperature );
        /* get random number seed from time if available 
            otherwise ask for a seed */
        printf( "Type number of configurations to average over:\n");
        ns = scanf( "%d", &nwobble );
        if( nwobble < 1 ) nwobble = 1;
        ltime = (long) time( NULL );
        iseed = (unsigned) ltime;
        if( ltime == -1 ) {
            printf("Type initial seed for random number generator:\n");
            ns = scanf("%ld", &iseed);
        } else {
            printf( "Random number seed initialized to %ld\n", iseed );
        }
        printf( "Type source size (FWHM in Ang.):\n" );
        ns = scanf( "%lf", &sourceFWHM );
    } else {
        temperature = 0.0F;
        nwobble = 1;
        sourceFWHM = 0.0;
    }
    /* convert FWHM to standard deviation 
            by dividing by 2*sqrt(2*ln(2)) */
    sourcesize = sourceFWHM / 2.354820045;

    timer = cputim();   /* get initial CPU time */
#ifdef USE_OPENMP
    walltimer = walltim();  /* wall time for opneMP */
#endif

    /*  calculate relativistic factor and electron wavelength */
    wavlen = (float) wavelength( keV );
    printf("electron wavelength = %g Angstroms\n", wavlen);

    /*---- read in specimen coordinates and scattering factors----- */

    natom = ReadXYZcoord( filein, ncellx, ncelly, ncellz,
        &ax, &by, &cz, &Znum, &xa, &ya, &za, &occ, &wobble,
        description, NCMAX );

    printf("%d atomic coordinates read in\n", natom );
    printf("%s", description );

    printf("Lattice constant a,b,c = %12.4f, %12.4f, %12.4f\n", ax,by,cz);
    
#ifdef USE_OPENMP
    /*  force LUT init. to avoid redundant init in parallel form */ 
    rsq = 0.5;  /* arbitrary position */   
    for( i=0; i<natom; i++) vz =  vzatomLUT( Znum[i], rsq );
#endif

    /* calculate thickness levels to save (1D mode) or check range (2D mode) */
    if( l1d == 1 ) {
        if( lxzimage == 1 ) {
            /* save all thickness levels  */
            nThick = (int) ( cz/deltaz + 0.5 );
            ThickSave = (double*) malloc1D( nThick, sizeof(double), "ThickSave");
            for( it=0; it<nThick; it++) {
                ThickSave[it] = deltaz*(it+1);
            }
        } else {
            nThick = 1;
            ThickSave = (double*) malloc1D( nThick, sizeof(double), "ThickSave");
            ThickSave[0] = cz;
        }
        printf( "save up to %d thickness levels\n", nThick );  /* diagnostic */
    } else {
        ThickSave[nThick-1] = cz;  /*  always save the last level */
        for( it=0; it<(nThick-1); it++) 
        if( (ThickSave[it] < 0.0) || (ThickSave[it] > cz) ) {
            printf("Bad thickness level = %g A, allowed range= "
                "0.0 to %f A\n", ThickSave[it], cz );
            exit( 0 );
        }
    }  /* end if( l1d == ... */

    if( lwobble == 0 ) {
        printf( "Sorting atoms by depth...\n");
        sortByZ( xa, ya, za, occ, Znum, natom );
    }
    /* to add random offsets */
    xa2 = (float*) malloc1D( natom, sizeof(float),  "xa2" );
    ya2 = (float*) malloc1D( natom, sizeof(float), "ya2" );
    za2 = (float*) malloc1D( natom, sizeof(float), "za2" );
    Znum2 = (int*) malloc1D( natom, sizeof(int), "Znum2" );
    occ2 = (float*) malloc1D( natom, sizeof(float), "occ2" );

    /*  calculate the total specimen volume and echo */
    xmin = xmax = xa[0];
    ymin = ymax = ya[0];
    zmin = zmax = za[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( xa[i] < xmin ) xmin = xa[i];
        if( xa[i] > xmax ) xmax = xa[i];
        if( ya[i] < ymin ) ymin = ya[i];
        if( ya[i] > ymax ) ymax = ya[i];
        if( za[i] < zmin ) zmin = za[i];
        if( za[i] > zmax ) zmax = za[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }
    printf("Total specimen range is\n %g to %g in x\n"
           " %g to %g in y\n %g to %g in z\n", xmin, xmax,
           ymin, ymax, zmin, zmax );
    if( lwobble == 1 )
        printf("Range of thermal rms displacements (300K) = %g to %g\n",
            wmin, wmax );

    /*  check for valid scan coordinates  */

    if( (xi < 0.0) || (xi > ax) ||
        (xf < 0.0) || (xf > ax) ||
        (yi < 0.0) || (yi > by) ||
        (yf < 0.0) || (yf > by) ) {
            printf("WARNING: Coordinates out of range; will be made periodic.\n");
            printf("xi,xf,yi,yf= %f, %f, %f, %f\n", xi, xf, yi, yf );
    }

    /*  check that requested probe size is not bigger 
        than transmission function size (or too small)
    */
    if( (nxprobe > nx) || (nxprobe < 2) ) {
        nxprobe = nx;
        printf("Probe size reset to nx = %d\n", nxprobe);
    }

    if( (nyprobe > ny) || (nyprobe < 2) ) {
        nyprobe = ny;
        printf("probe size reset to ny = %d\n", nyprobe);
    }

    /*  calculate spatial frequencies for future use
        (one set for transmission function and one for probe
        wavefunction)
    NOTE: zero freg is in the bottom left corner and
        expands into all other corners - not in the center
        this is required for FFT - don't waste time rearranging

    remember : the x values are the same for both sets
    
    x2, y2 are used for confocal
    
    */

    kx  = (float*) malloc1D( nx, sizeof(float), "kx" );
    ky  = (float*) malloc1D( ny, sizeof(float), "ky" );
    kx2 = (float*) malloc1D( nx, sizeof(float), "kx2" );
    ky2 = (float*) malloc1D( ny, sizeof(float), "ky2" );
    xp  = (float*) malloc1D( nx, sizeof(float), "x2" );
    yp  = (float*) malloc1D( ny, sizeof(float), "y2" );

    freqn( kx, kx2, xp, nx, ax );
    freqn( ky, ky2, yp, ny, by );

    kxp  = (float*) malloc1D( nxprobe, sizeof(float), "kxp" );
    kyp  = (float*) malloc1D( nyprobe, sizeof(float), "kyp" );
    kxp2 = (float*) malloc1D( nxprobe, sizeof(float), "kxp2" );
    kyp2 = (float*) malloc1D( nyprobe, sizeof(float), "kyp2" );
    
    freqn( kxp, kxp2, xp, nxprobe, dxp = ax*((double)nxprobe)/nx );
    freqn( kyp, kyp2, yp, nyprobe, dyp = by*((double)nyprobe)/ny );

    /* impose anti-aliasing bandwidth limit on transmission functions */

    sum = ((double)nx)/(2.0*ax);
    k2maxp = ((double)ny)/(2.0*by);
    if( sum < k2maxp ) k2maxp = sum;
    k2maxp= BW * k2maxp;
    printf("Bandwidth limited to a real space resolution of %f Angstroms\n",
                     1.0F/k2maxp);
    printf("   (= %.2f mrad)  for symmetrical anti-aliasing.\n",
         wavlen*k2maxp*1000.0F);
    k2maxp = k2maxp * k2maxp;

    /*  allocate some more arrays and initialize propagator */

    propxr = (float*) malloc1D( nxprobe, sizeof(float), "propxr" );
    propxi = (float*) malloc1D( nxprobe, sizeof(float), "propxi" );
    propyr = (float*) malloc1D( nyprobe, sizeof(float), "propyr" );
    propyi = (float*) malloc1D( nyprobe, sizeof(float), "propyi" );

    /* calculate propagator functions with probe sample size
        impose anti-aliasing bandwidth limit  */
    tctx = 2.0 * tan(ctiltx);
    tcty = 2.0 * tan(ctilty);

    scale = pi * deltaz;
    for( ix=0; ix<nxprobe; ix++) {
        w = scale * ( kxp2[ix] * wavlen - kxp[ix]*tctx );
        propxr[ix]= (float)  cos(w);
        propxi[ix]= (float) -sin(w);
    }

    for( iy=0; iy<nyprobe; iy++) {
        w = scale * ( kyp2[iy] * wavlen - kyp[iy]*tcty );
        propyr[iy]= (float)  cos(w);
        propyi[iy]= (float) -sin(w);
    }

    /*   calculate number of pixels in the probe and obj. apert. */
    k2maxa = apert1 * 0.001/wavlen;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 * 0.001/wavlen;
    k2maxb = k2maxb * k2maxb;
    
    nbeamp = nbeampo = 0;
    for( iy=0; iy<nyprobe; iy++)
    for( ix=0; ix<nxprobe; ix++) {
        k2 = kyp2[iy] + kxp2[ix];
        if( k2 < k2maxp ) nbeamp++;
        if( (k2 >= k2maxa) && (k2 <= k2maxb) ) nbeampo++;
    }

    printf("Number of symmetrical anti-aliasing "
           "beams in probe = %ld\n", nbeamp);
    printf("Number of beams in probe aperture = %ld\n", nbeampo);
    if( nbeamp < 200 ) {
        printf("WARNING: the probe is under sampled, this is a bad idea...\n");
    }
    if( nbeampo < 100 ) {
        printf("WARNING: the probe aperture is under sampled, this is a bad idea...\n");
        exit( 0 );
    }

    /*  convert aperture dimensions */

    k2min = (double*) malloc1D( ndetect, sizeof(double), "k2min" );
    k2max = (double*) malloc1D( ndetect, sizeof(double), "k2max" );

    for( idetect=0; idetect<ndetect; idetect++) {
        if( ADF == collectorMode[idetect] ) {
            k2max[idetect] = 0.001 * almax[idetect]/wavlen;
            k2max[idetect] = k2max[idetect] * k2max[idetect];
            k2min[idetect] = 0.001 * almin[idetect]/wavlen;
            k2min[idetect] = k2min[idetect] * k2min[idetect];
        } else if( CONFOCAL == collectorMode[idetect] ) {
            k2max[idetect] = almax[idetect] * almax[idetect];
            k2min[idetect] = almin[idetect] * almin[idetect];
        }
    }

    /*  init the min/max record of total integrated intensity */

    totmin =  10.0;
    totmax = -10.0;
    detect  = (double***) malloc3D( nThick, ndetect, nprobes,
        sizeof(double), "detect" );
    detect2 = (double***) malloc3D( nThick, ndetect, nprobes,
        sizeof(double), "detect2" );
    sums = (double*) malloc1D( nprobes, sizeof(double), "sums" ); 
    rmin = (float**) malloc2D( nThick, ndetect, sizeof(float), "rmin" );
    rmax = (float**) malloc2D( nThick, ndetect, sizeof(float), "rmax" );

    /* allocate probe wave function and transmission function arrays
         using FFTW 3 data structure and make FFTW plan */

    probe = (fftwf_complex**) malloc( nprobes * sizeof( fftwf_complex* ) );
        if( NULL == probe ) {
        printf("Cannot allocate probe array\n");
        exit( EXIT_FAILURE );
    }
    for( ip=0; ip<nprobes; ip++){
        probe[ip] = (fftwf_complex*) fftwf_malloc( nxprobe*nyprobe * sizeof(fftwf_complex) );
        if( NULL == probe[ip] ) {
            printf("Cannot allocate probe array\n");
            exit( EXIT_FAILURE );
        }
    }
    /* remember FFTW uses opposite sign convention */
    planPi = fftwf_plan_dft_2d( nxprobe, nyprobe, probe[0], probe[0], 
         FFTW_FORWARD, FFTW_MEASURE );  /* inverse in place */
    planPf = fftwf_plan_dft_2d( nxprobe, nyprobe, probe[0], probe[0], 
         FFTW_BACKWARD, FFTW_MEASURE );  /* forward in place */

    trans = (fftwf_complex*) fftwf_malloc( nx*ny * sizeof(fftwf_complex) );
    if( NULL == trans ) {
        printf("Cannot allocate trans array\n");
        exit( EXIT_FAILURE );
    }
    /* remember FFTW uses opposite sign convention */
    planTi = fftwf_plan_dft_2d( nx, ny, trans, trans, 
        FFTW_FORWARD, FFTW_MEASURE );  /* inverse in place */
    planTf = fftwf_plan_dft_2d( nx, ny, trans, trans, 
        FFTW_BACKWARD, FFTW_MEASURE );  /* forward in place */

    if( lpacbed == TRUE ) {
        pacbedPix = (float**) malloc2D( nxprobe, nyprobe, sizeof(float), "pacbedPix" );
        for( ix=0; ix<nxprobe; ix++) for( iy=0; iy<nyprobe; iy++)
                pacbedPix[ix][iy] = 0;
    }

/* ------------- start here for a full image output -------------- */
/*
  do one whole line at once NOT the whole image (which may be huge)
*/
    if( l1d == 0 ) {
       printf("output file size in pixels is %d x %d\n",
          nxout, nyout );
       if( nprobes != nyout ) {
           printf( "Error, nprobes=%d must be the same as"
               "nyout=%d, in image mode.\n", nprobes, nyout );
           exit( 0 );
       }
       /* double up first index to mimic a 4D array */
       pixr = (float***) malloc3D( ndetect*nThick, nxout, nyout,
           sizeof(float), "pixr"  );
       for( i=0; i<(nThick*ndetect); i++) {
           for( ix=0; ix<nxout; ix++)
           for( iy=0; iy<nyout; iy++)
            pixr[i][ix][iy] = 0.0F;
       }

       /*  iterate the multislice algorithm proper for each
           position of the focused probe */

       dx = (xf-xi)/((double)(nxout-1));
       dy = (yf-yi)/((double)(nyout-1));
       x = (double*) malloc1D( nprobes, sizeof(double), "x" );
       y = (double*) malloc1D( nprobes, sizeof(double), "y" );

        /*  add random thermal displacements 
               scaled by temperature if requested 
            remember that initial wobble is at 300K for
               each direction */

        for( iwobble=0; iwobble<nwobble; iwobble++) {
            if( lwobble == 1 ){
                scale = (float) sqrt(temperature/300.0) ;
                for( i=0; i<natom; i++) {
                    xa2[i] = xa[i] + 
                        (float)(wobble[i]*rangauss(&iseed)*scale);
                    ya2[i] = ya[i] + 
                        (float)(wobble[i]*rangauss(&iseed)*scale);
                    za2[i] = za[i] + 
                            (float)(wobble[i]*rangauss(&iseed)*scale);
                    occ2[i] = occ[i];
                    Znum2[i] = Znum[i];
                }
                sortByZ( xa2, ya2, za2, occ2, Znum2, natom );
                printf("configuration # %d\n", iwobble+1 );
                printf( "The new range of z is %g to %g\n",
                    za2[0], za2[natom-1] );
            } else for( i=0; i<natom; i++) {
                xa2[i] = xa[i];
                ya2[i] = ya[i];
                za2[i] = za[i];
                occ2[i] = occ[i];
                Znum2[i] = Znum[i];
            }
            zmin = za2[0];  /* reset zmin/max after wobble */
            zmax = za2[natom-1];
    
            for( ix=0; ix<nxout; ix++) {
    
                for( iy=0; iy<nyout; iy++) {
                    x[iy] = xi + dx * ((double) ix)
                        + sourcesize * rangauss(&iseed);
                    y[iy] = yi + dy * ((double) iy)
                        + sourcesize * rangauss(&iseed);
                    x[iy] = periodic( x[iy], ax );   /* put back in supercell */
                    y[iy] = periodic( y[iy], by );   /* if necessary */
                }

                STEMsignals( x, y, nyout, param, multiMode, detect, ndetect, 
                    ThickSave, nThick, sums );
                for( iy=0; iy<nyout; iy++) {
                    if( sums[iy] < totmin ) totmin = sums[iy];
                    if( sums[iy] > totmax ) totmax = sums[iy];
                    for( it=0; it<nThick; it++){
                        for( idetect=0; idetect<ndetect; idetect++)
                        pixr[idetect + it*ndetect][ix][iy] += (float)
                            (detect[it][idetect][iy]/((double)nwobble));
                    }
                    if( sums[iy] < 0.9) 
                        printf("Warning integrated intensity too small, = "
                        "%g at x,y= %g, %g\n", sums[iy], x[iy], y[iy] );
                    if( sums[iy] > 1.1) 
                        printf("Warning integrated intensity too large, = "
                        "%g at x,y= %g, %g\n", sums[iy], x[iy], y[iy] );
                }

                /*   sum position averaged CBED if requested 
                     - assume probe still left from stemsignal()  */
                if( lpacbed == TRUE ) {
                    for( iy=0; iy<nyout; iy++) {
                        for( ix2=0; ix2<nxprobe; ix2++)
                        for( iy2=0; iy2<nyprobe; iy2++) {
                            prr = probe[iy][iy2 + ix2*nyprobe][0];
                            pri = probe[iy][iy2 + ix2*nyprobe][1];
                            pacbedPix[ix2][iy2] += (prr*prr + pri*pri);
                       }
                    }
                }   /*  end if( lpacbed.... */
            
            } /* end for(ix...) */
    
        } /* end for(iwobble... ) */

        /*  output data files  */
        for( it=0; it<nThick; it++)
        for( i=0; i<ndetect; i++) {
            rmin[it][i] = rmax[it][i] = pixr[i+it*ndetect][0][0];
            for( ix=0; ix<nxout; ix++)
            for( iy=0; iy<nyout; iy++) {
                temp = pixr[i+it*ndetect][ix][iy];
                if( temp < rmin[it][i] )rmin[it][i] = (float) temp;
                if( temp > rmax[it][i] )rmax[it][i] = (float) temp;
            }
        }

        /*  directory file listing parameters for each image file */
        sprintf( fileout, "%s.txt", fileoutpre );
        fp = fopen( fileout, "w+" );
        if( fp == NULL ) {
            printf("Cannot open output file %s.\n", fileout );
            exit( 0 );
        }
     
       fprintf(fp, "C\n");
       fprintf(fp, "C   output of autostem version %s\n", version);
       fprintf(fp, "C\n");
       fprintf(fp, "C   nslice= %d\n", nslice);
       fprintf(fp,"deltaz= %g, file in= %s\n", deltaz, filein );
       fprintf(fp, "V0= %g, Cs3= %g, Cs5= %g, df= %g\n", keV, Cs3, Cs5, df );
       fprintf(fp, "Apert= %g mrad to %g mrad\n", apert1, apert2 );
       if( doConfocal == TRUE ) {
          fprintf(fp, "confocal lens Cs3= %g, Cs5= %g, df= %g\n",  Cs3C, Cs5C, dfC );
          fprintf(fp, "confocal lens apert= %g mrad to %g mrad\n", apert1C, apert2C );
          fprintf(fp, "confocal dfa2C= %g, dfa2phiC= %g, dfa3C= %g, dfa3phiC=%g\n",
                   dfa2C, dfa2phiC, dfa3C, dfa3phiC );
       }
       fprintf(fp, "Crystal tilt x,y= %lg, %lg\n", ctiltx,ctilty);

       for(  idetect=0; idetect<ndetect; idetect++) {
            if( ADF == collectorMode[idetect]  )
                fprintf(fp, "Detector %d, Almin= %g mrad, Almax= %g mrad\n",
                            idetect, almin[idetect], almax[idetect] );
            else if( CONFOCAL == collectorMode[idetect] )
                fprintf(fp, "Detector %d, cmin= %g Angst, cmax= %g Angst.\n",
                            idetect, almin[idetect], almax[idetect] );
       }

       fprintf(fp, "ax= %g A, by= %g A, cz= %g A\n", ax,by,cz);
       fprintf(fp, "Number of symmetrical anti-aliasing "
           "beams in probe wave function= %ld\n", nbeamp );
       fprintf(fp, "with a resolution (in Angstroms) = %g\n",
           1.0/sqrt(k2maxp) );
       if( lwobble == 1 ) {
          fprintf( fp,
           "Number of thermal configurations = %d\n", nwobble );
              fprintf( fp, "Source size = %g Ang. (FWHM) \n", sourceFWHM );
       }
       fprintf( fp, "\n" );

        /*  store params plus min and max */
        param[pIMAX]    = 0.0F;
        param[pIMIN]    = 0.0F;
        param[pXCTILT]  = (float) ctiltx;
        param[pYCTILT]  = (float) ctilty;
        param[pDEFOCUS] = (float) df;
        param[pDX]      = (float) dx;
        param[pDY]      = (float) dy;
        param[pENERGY]  = (float) keV;
        param[pOAPERT]  = (float) apert2;
        param[pWAVEL]   = (float) wavlen;
        param[pNSLICES] = (float) -1.0; /* ??? ejk 19-jul-1998 */

        for( i=0; i<NPARAM; i++) myFile.setParam( i, param[i] );
        myFile.setnpix( 1 );
        myFile.resize( nxout, nyout );
        aimin = aimax = 0.0F;

        for( it=0; it<nThick; it++)
        for( i=0; i<ndetect; i++) {
            sprintf( fileout, "%s%03d%03d.tif", fileoutpre, i, it );
            printf("%s: output pix range : %g to %g\n", fileout, 
                rmin[it][i], rmax[it][i]);
            myFile.setParam( pRMAX, rmax[it][i] );
            myFile.setParam( pRMIN, rmin[it][i] );
            myFile.setParam( pMINDET, (float) ( almin[i] * 0.001 ) );
            myFile.setParam( pMAXDET, (float) ( almax[i] * 0.001 ) );
            for( ix=0; ix<nxout; ix++) for( iy=0; iy<nyout; iy++)  
                myFile( ix, iy ) = pixr[i+it*ndetect][ix][iy];
            if( myFile.write( fileout, rmin[it][i], rmax[it][i], aimin, aimax,
                (float) dx, (float) dy ) != 1 ) {
                    printf("Cannot write output file %s.\n", fileout );
            }

            if( ADF == collectorMode[i]  )
                fprintf(fp,"file: %s, detector= %g to %g mrad, "
                        "thicknes= %g A, range= %g to %g\n", fileout,
                        almin[i], almax[i], ThickSave[it], rmin[it][i], rmax[it][i]);
            else if( CONFOCAL == collectorMode[i] )
                fprintf(fp,"file: %s, detector= %g to %g Angst., "
                        "thicknes= %g A, range= %g to %g\n", fileout,
                        almin[i], almax[i], ThickSave[it], rmin[it][i], rmax[it][i]);
        }  /*  end for(i=... */

        /*   save pos. aver. CBED if needed */
        if( lpacbed == TRUE ) {
            myFile.setnpix( 1 );
            invert2D( pacbedPix, nxprobe, nyprobe );  /*  put zero in middle */
            nx1 =    nxprobe / 6;
            nx2 = (5*nxprobe) / 6;  /*  cut out center portion without anti-aliasing zeros */
            ny1 =    nyprobe / 6;
            ny2 = (5*nyprobe) / 6;
            nxout2 = nx2 - nx1 + 1;
            nyout2 = ny2 - ny1 + 1;
	    if( (nxout2<1) || (nyout2<1) ) {
		    nx1 = ny1 = 0;
		    nx2 = nxout2 = nxprobe;
		    ny2 = nyout2 = nyprobe;
	    }
            myFile.resize( nxout2, nyout2 );
            aimin = aimax = 0.0F;
	    scalef = (float) (1.0/(nxout2*nyout2) );
            ixo = 0;
            for( ix2=nx1; ix2<=nx2; ix2++) {
                iyo = 0;
                for( iy2=ny1; iy2<=ny2; iy2++) {
                    myFile(ixo,iyo++) = scalef * pacbedPix[ix2][iy2];
                }  ixo++;
            }
            rmin0 = myFile.min(0);
            rmax0 = myFile.max(0);
            dxp = 1.0/dxp;
            dyp = 1.0/dyp;
            myFile.setParam( pRMAX, rmax0 );
            myFile.setParam( pRMIN, rmin0 );
            myFile.setParam( pDX, (float) dxp);
            myFile.setParam( pDY, (float) dyp );
            printf( "pos. averg. CBED (unaliased) size %d x %d pixels\n"
                " and range (arb. units): %g to %g\n", nxout2, nyout2, rmin0, rmax0 );
            if( myFile.write( pacbedFile, rmin0, rmax0, aimin, aimax,
                (float) dxp, (float) dyp ) != 1 ) {
                printf("Cannot write output file %s.\n", pacbedFile );
            }
        }   /*  end if( lpacbed.... */

        fclose( fp );

/* ------------- start here for 1d line scan ---------------- */

    } else if ( l1d == 1 ) {

       if( lpacbed == TRUE ) {
           printf("warning: cannot do pos. aver. CBED in 1d\n");
       }
       dx = (xf-xi)/((double)(nout-1));
       dy = (yf-yi)/((double)(nout-1));
       x = (double*) malloc1D( nprobes, sizeof(double), "x" );
       y = (double*) malloc1D( nprobes, sizeof(double), "y" );
            for( ip=0; ip<nout; ip++) {
                for( it=0; it<nThick; it++)
                for( idetect=0; idetect<ndetect; idetect++)
                    detect[it][idetect][ip] = 0.0;
       }

       /*  add random thermal displacements scaled by temperature
            if requested 
        remember that initial wobble is at 300K for each direction */

       for( iwobble=0; iwobble<nwobble; iwobble++) {
    
            if( lwobble == 1 ){
                scale = (float) sqrt(temperature/300.0) ;
                for( i=0; i<natom; i++) {
                    xa2[i] = xa[i] + 
                        (float)(wobble[i]*rangauss(&iseed)*scale);
                    ya2[i] = ya[i] + 
                        (float)(wobble[i]*rangauss(&iseed)*scale);
                    za2[i] = za[i] + 
                        (float)(wobble[i]*rangauss(&iseed)*scale);
                    occ2[i] = occ[i];
                    Znum2[i] = Znum[i];
                }
                printf("configuration # %d\n", iwobble+1 );
                sortByZ( xa2, ya2, za2, occ2, Znum2, natom );
                printf( "The new range of z is %g to %g\n",
                    za2[0], za2[natom-1] );
            } else for( i=0; i<natom; i++) {
                xa2[i] = xa[i];
                ya2[i] = ya[i];
                za2[i] = za[i];
                occ2[i] = occ[i];
                Znum2[i] = Znum[i];
            }
            zmin = za2[0];      /* reset zmin/max after wobble */
            zmax = za2[natom-1];
            for( ip=0; ip<nout; ip++) {
                x[ip] = xi + dx * ((double)ip)
                            + sourcesize * rangauss(&iseed);
                y[ip] = yi + dy * ((double)ip)
                            + sourcesize * rangauss(&iseed);
                x[ip] = periodic( x[ip], ax );   /* put back in supercell */
                y[ip] = periodic( y[ip], by );  /* if necessary */
            }
           
            STEMsignals( x, y, nprobes, param, multiMode, detect2, ndetect, 
                ThickSave, nThick, sums );
            for( ip=0; ip<nprobes; ip++) {
                if( sums[ip] < totmin ) totmin = sums[ip];
                if( sums[ip] > totmax ) totmax = sums[ip];
                for( it=0; it<nThick; it++){
                for( idetect=0; idetect<ndetect; idetect++)
                   detect[it][idetect][ip] += 
                        detect2[it][idetect][ip]/((double)nwobble);
                }
                if( sums[ip] < 0.9) 
                printf("Warning integrated intensity too small, = %g"
                    " at x,y= %g, %g\n", sums[ip], x[ip], y[ip] );
                if( sums[ip] > 1.1) 
                printf("Warning integrated intensity too large, = %g"
                    " at x,y= %g, %g\n", sums[ip], x[ip], y[ip] );
            }

       }  /* end for(iwobble... */

       /* ------ first output text data ---------------- */
       sprintf( fileout, "%s.dat", fileoutpre );
       printf("output file= %s\n", fileout);

       fp = fopen( fileout, "w+" );
       if( fp == NULL ) {
           printf("Cannot open output file %s.\n", fileout );
             exit( 0 );
       }

       fprintf(fp, "C\n");
       fprintf(fp, "C   output of autostem version %s\n", version);
       fprintf(fp, "C\n");
       fprintf(fp, "C   nslice= %d\n", nslice);
       fprintf(fp,"deltaz= %g, file in= %s\n", deltaz, filein );
       fprintf(fp, "V0= %g, Cs3= %g, Cs5= %g, df= %g\n", keV, Cs3, Cs5, df );
       fprintf(fp, "Apert= %g mrad to %g mrad\n", apert1, apert2 );
       if( doConfocal == TRUE ) {
          fprintf(fp, "confocal lens Cs3= %g, Cs5= %g, df= %g\n",  Cs3C, Cs5C, dfC );
          fprintf(fp, "confocal lens apert= %g mrad to %g mrad\n", apert1C, apert2C );
          fprintf(fp, "confocal dfa2C= %g, dfa2phiC= %g, dfa3C= %g, dfa3phiC=%g\n",
                   dfa2C, dfa2phiC, dfa3C, dfa3phiC );
       }
       fprintf(fp, "Crystal tilt x,y= %lg, %lg\n", ctiltx,ctilty);

       for(  idetect=0; idetect<ndetect; idetect++) {
            if( ADF == collectorMode[idetect]  )
                fprintf(fp, "Detector %d, Almin= %g mrad, Almax= %g mrad\n",
                            idetect, almin[idetect], almax[idetect] );
            else if( CONFOCAL == collectorMode[idetect] )
                fprintf(fp, "Detector %d, cmin= %g Angst, cmax= %g Angst.\n",
                            idetect, almin[idetect], almax[idetect] );
       }

       fprintf(fp, "ax= %g A, by= %g A, cz= %g A\n", ax,by,cz);
       fprintf(fp, "Number of symmetrical anti-aliasing "
           "beams in probe wave function= %ld\n", nbeamp );
       fprintf(fp, "with a resolution (in Angstroms) = %g\n",
           1.0/sqrt(k2maxp) );
       if( lwobble == 1 ) {
              fprintf( fp,
               "Number of thermal configurations = %d\n", nwobble );
                  fprintf( fp, "Source size = %g Ang. (FWHM) \n", sourceFWHM );
       }
       fprintf(fp, "C     x      y     signal\n");

       for( ip=0; ip<nprobes; ip++) {
           /* recalculate mean x,y without source size wobble */
           x[ip] = xi + dx * ((double)ip);
           y[ip] = yi + dy * ((double)ip);
           fprintf(fp, "%14.7g %14.7g", x[ip], y[ip]);
           for(i=0; i<ndetect; i++) 
               fprintf(fp, "%14.7g", detect[nThick-1][i][ip] );
           fprintf(fp, "\n");
       }

       fclose( fp );

       /* ------ next output xz image data ---------------- */

       if( lxzimage == 1 ) {
    
            /*  directory file listing parameters for each image file */
            sprintf( fileout, "%s.txt", fileoutpre );
            fp = fopen( fileout, "w+" );
            if( fp == NULL ) {
                printf("Cannot open output file %s.\n", fileout );
                exit( 0 );
            }
    
            /*  store params plus min and max */
            param[pIMAX]    = 0.0F;
            param[pIMIN]    = 0.0F;
            param[pXCTILT]  = (float) ctiltx;
            param[pYCTILT]  = (float) ctilty;
            param[pDEFOCUS] = (float) df;
            param[pDX]  = (float) dx;
            param[pDY]  = (float) dy;
            param[pENERGY]  = (float) keV;
            param[pOAPERT]  = (float) apert2;
            param[pWAVEL]   = (float) wavlen;
            param[pNSLICES] = (float) -1.0; /* ??? ejk 19-jul-1998 */
    
            myFile.setnpix( 1 );
            myFile.resize( nprobes, nThick );
            aimin = aimax = 0.0F;
            for( i=0; i<NPARAM; i++) myFile.setParam( i, param[i] );
  
            for( idetect=0; idetect<ndetect; idetect++){
                sprintf( fileout, "%s%03d.tif", fileoutpre, idetect );
                printf("output file= %s\n", fileout);
    
                /* convert to float and fix pixel order */
                pmin = pmax = (float) detect[0][idetect][0];
                for( ix=0; ix<nprobes; ix++)
                for( iy=0; iy<nThick; iy++) {
                    temp = myFile( ix, iy ) = (float) detect[iy][idetect][ix];
                    if( temp < pmin )pmin = temp;
                    if( temp > pmax )pmax = temp;
                }
    
                printf("%s: output pix range : %g to %g\n", fileout, pmin, pmax);
                myFile.setParam(pRMAX, pmax );
                myFile.setParam(pRMIN, pmin );
                if( collectorMode[idetect] == ADF ) {
                    myFile.setParam(pMINDET, (float) ( almin[idetect] * 0.001 ) );
                    myFile.setParam(pMAXDET, (float) ( almax[idetect] * 0.001 ) );
                } else if( collectorMode[idetect] == CONFOCAL ) {
                    myFile.setParam(pMINDET, (float) almin[idetect] );
                    myFile.setParam(pMAXDET, (float) almax[idetect] );
                }
 
                if(  myFile.write( fileout, pmin, pmax,
                    aimin, aimax, (float) dx, (float) dy ) != 1 ) {
                        printf("Cannot write output file %s.\n", fileout );
                }
                fprintf(fp,"file: %s, detector= %g to %g mrad, range= %g to %g\n",
                    fileout, almin[idetect], almax[idetect], pmin, pmax);
            }
            fclose( fp );

       }  /* end if( lxzimage==1... */

    } /* end if( l1d.. ) */

    printf("Number of symmetrical anti-aliasing "
           "beams in trans. function = %ld\n", nbeamt);

    /*  echo min/max of total integrated intensity */
    printf("The total integrated intensity range was:\n");
    printf("   %g to %g\n\n",  totmin, totmax );

    printf("CPU time = %g sec.\n", cputim()-timer);
#ifdef USE_OPENMP
    printf("wall time = %g sec.\n", walltim() - walltimer);
#endif

    return( 0 );

}  /* end main() */

/*------------------------ periodic() ---------------------*/
/*
     make probe positions periodic in the supercell
     in case some wobble off the edge with source size of user excess

    pos = input position (x or y);
    size = supercell size ( 0 to size)

    return positive value  0 <= x < size
*/
double periodic( double pos, double size )
{
    double x=pos;
    while( x < 0 ) x += size;
    x = fmod( x, size );
    return( x );
}

/*------------------------ STEMsignals() ---------------------*/
/*

  NOTE: this is NOT the same as STEMsignal() in stemslice.c

  subroutine to calculate the stem signal arising from a given
  probe position

  iterate the multislice algorithm proper for each position of
  the focused probe

  This version uses massive amounts of memory to avoid
  recalculating the transmission functions more than necessary

   note zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for fft - don't waste time rearranging

     change to propagate thru whole unit cell not just
     range of atoms on 14-sep-2006 ejk
     add multipole aberrations 9-may-2011 ejk

  x[],y[]     = real positions of the incident probe
  npos        = int number of positions
  param[]     = parameters of probe
  multiMode   = flag to add multipole aberrations
  detect[][][]= real array to get signal into each detector
            for each probe position and thickness
  ndetect     = number of detector geometries
  ThickSave[] = thicknesses at which to save data (other than the last)
  nThick      = number of thickness levels (including the last)
  sum         = real total integrated intensity
  
  the assumed global variables are:
  
  nxprobe,nyprobe   = int size of probe wavefunction in pixels
  nx,ny         = int size of transmission function in pixels
  layer[]       = int array with slice layer indecies
  prober[][], probei[][] = float real,image probe wavefunction
  transr[][], transi[][] = float real,imag transmission function
  propxr[][], propxi[][] = float real,imag propagator vs x
  propyr[][], propyi[][] = float real,imag propagator vs y
  ax,by,cz      = float unit cell size in Angs
  kxp[], kyp[]      = float spatial frequencies vs x, y
  kxp2[], kyp2[]    = float square of kxp[], kyp[]
  xp[], yp[]        = float real space positions in probe (confocal)
  apert1, apert2    = double min,max objective aperture (mrad)
  k2maxp            = double max spatial freq of probe squared
  pi                = double constant PI
  wavlen            = double electron wavelength in Angs
  df                = double defocus (in Ang)
  Cs3,Cs5           = double spherical aberration (in mm)

  xa[],ya[],za[]    = atom coordinates
  occ[]         = atomic occupancy
  Znum[]        = atomic numbers
  natom         = number of atoms
  deltaz        = slice thickness
  v0            = beam energy
  nbeamt        = number of beams in transmission function
  zmin, zmax        = range of z coord. of the atoms
  nslice        = number of slices
  
    NOTE:  too many thing come in as globals, but...

*/
void STEMsignals( double x[], double y[], int npos, float p[],
         int multiMode, double ***detect, int ndetect,
         double ThickSave[], int nThick, double sum[] )
{
    int ix, iy, ixt, iyt, idetect, *ixoff, *iyoff, ixmid, iymid;
    int istart, na, ip, i, j, jt, it;

    long nxprobel, nyprobel, nxl, nyl;

    float scale, prr, pri, tr, ti;
    fftwf_complex *cpix;            /* complex confocal image */

    double *xoff, *yoff, chi0, chi1, k2maxa, k2maxb,
        w, k2, phi, phirms, alx, aly;
    double sum0, sum1, delta, zslice, totalz;
    
    /* extra for confocal */
    float hr, hi;
    double chi2C, chi3C, k2maxaC, k2maxbC, r2, rx2;

    /* ------ make sure x,y are ok ------ */

    for( ip=0; ip<npos; ip++) {
        if( (x[ip] < 0.0) || (x[ip] > ax) ||
            (y[ip] < 0.0) || (y[ip] > by) ) {
            sum[ip] = -1.2345;
            printf("bad x=%f,y=%f in STEMsignals()\n", x[ip], y[ip]);
            return;
        }
    }

    /*  generate a probe at position x,y which must be inside the
        0->ax and 0->by region
        normalize to have unit integrated intensity
        probe position (x,y) = (0,0) = lower left corner

    NOTE: The probe wavefunction is nxprobe x nyprobe which may
        be smaller than the transmission function.
        Position the probe near the center of the nxprobe x nyprobe
        region because it does not wrap around properly
        (i.e. manually wrap around in the nx x ny region if
         the probe is near a boundary of the nx x ny region))
    */

    ixmid = nxprobe/2;
    iymid = nyprobe/2;
    chi1 = pi * wavlen;
    k2maxa = apert1 * 0.001/wavlen;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 * 0.001/wavlen;
    k2maxb = k2maxb * k2maxb;
    
    /* extra for confocal */
    chi2C = 0.5 * Cs3C * 1.e7*wavlen*wavlen;
    chi3C = Cs5C * 1.0e7 * wavlen*wavlen*wavlen*wavlen /3.0;
    k2maxaC = apert1C * 0.001/wavlen;
    k2maxaC = k2maxaC *k2maxaC;
    k2maxbC = apert2C * 0.001/wavlen;
    k2maxbC = k2maxbC * k2maxbC;

    ixoff = (int*) malloc1D( npos, sizeof(int), "ixoff" );
    iyoff = (int*) malloc1D( npos, sizeof(int), "iyoff" );
    xoff = (double*) malloc1D( npos, sizeof(double), "xoff" );
    yoff = (double*) malloc1D( npos, sizeof(double), "yoff" );

    /* ------- calculate all of the probe wave functions at once ------
        to reuse the transmission functions which takes a long
        time to calculate*/

/*  paralleling this loop has little effect */
#pragma omp parallel for private(ix,iy,j,sum0,k2,w,chi0,scale,tr,ti,alx,aly) 
    for( ip=0; ip<npos; ip++) {
        ixoff[ip] = (int) floor( x[ip]*((double)nx) / ax ) - ixmid;
        xoff[ip]  = x[ip] - ax*((double)ixoff[ip])/((double)nx);

        iyoff[ip] = (int) floor( y[ip]*((double)ny) / by ) - iymid;
        yoff[ip]  = y[ip] - by*((double)iyoff[ip])/((double)ny);

        sum0 = 0.0;
        for( ix=0; ix<nxprobe; ix++) {
            j = ix*nyprobe;
            alx = wavlen * kx[ix];  /* x component of angle alpha */
            for( iy=0; iy<nyprobe; iy++) {
                aly = wavlen * ky[iy];  /* y component of angle alpha */
                k2 = kxp2[ix] + kyp2[iy];
                if( (k2 >= k2maxa) && (k2 <= k2maxb) ) {
                    w = 2.*pi* ( xoff[ip]*kxp[ix] + yoff[ip]*kyp[iy] );
                    chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode );
                    chi0= - chi0 + w;
                    probe[ip][iy + j][0] = tr = (float) cos( chi0 );
                    probe[ip][iy + j][1] = ti = (float) sin( chi0 );
                    sum0 += (double) (tr*tr + ti*ti);
                } else {
                    probe[ip][iy + j][0] = 0.0F;
                    probe[ip][iy + j][1] = 0.0F;
                }
            }
        }  /* end for( ix... */

        scale = (float) ( 1.0/sqrt(sum0) );
        for( ix=0; ix<nxprobe; ix++) {
            j = ix*nyprobe;
            for( iy=0; iy<nyprobe; iy++) {
                probe[ip][j  ][0] *= scale;
                probe[ip][j++][1] *= scale;
            }
        }

    }  /* end for( ip...) */

    /* -------- transmit thru nslice layers ------------------------
        beware ixoff,iyoff must be with one nx,ny
            of the array bounds
    */

    nxprobel = (long) nxprobe;
    nyprobel = (long) nyprobe;

    nxl = (long) nx;
    nyl = (long) ny;
    
    scale = 1.0F / ( ((float)nx) * ((float)ny) );

    zslice = 0.75*deltaz;  /*  start a little before top of unit cell */
    istart = 0;
    nslice = 0;
    it = 0;     /* thickness level index */

    if( zmax > cz ) totalz = zmax;
        else totalz = cz;
    printf( "specimen range is 0 to %g Ang.\n", totalz );
    
    /* range of unit cell */
    while(  (zslice < (totalz+0.25*deltaz)) || (istart<natom) ) {

       /* find range of atoms for current slice */
       na = 0;
       for(i=istart; i<natom; i++) 
            if( za2[i] < zslice ) na++; else break;

       printf( "slice ending at z= %g Ang. with %d atoms\n", zslice, na );

       /* calculate transmission function and bandwidth limit */
       if( na > 0 ) trlayer( &xa2[istart], &ya2[istart], &occ2[istart],
            &Znum2[istart], na, (float)ax, (float)by, (float)keV,
            trans, nxl, nyl, &phirms, &nbeamt, (float) k2maxp );

            /*----- one multislice trans/prop cycle for all probes ---- */
#pragma omp parallel for private(ix,iy,ixt,iyt,j,jt,prr,pri)
       for( ip=0; ip<npos; ip++) {
           /* apply transmission function if there are atoms in this slice */
           if( na > 0 ) {
                fftwf_execute_dft( planPi, probe[ip], probe[ip] );
                scaleW( probe[ip], nxprobe, nyprobe );
    
                for( ix=0; ix<nxprobe; ix++) {
                    ixt = ix + ixoff[ip];
                    if( ixt >= nx ) ixt = ixt - nx;
                    else if( ixt < 0 ) ixt = ixt + nx;
                    j = ix*nyprobe;
                    jt = ixt*ny;
                    for( iy=0; iy<nyprobe; iy++) {
                        iyt = iy + iyoff[ip];
                        if( iyt >= ny ) iyt = iyt - ny;
                        else if( iyt < 0 ) iyt = iyt + ny;
                        prr = probe[ip][iy + ix*nyprobe][0];
                        pri = probe[ip][iy + ix*nyprobe][1];
                        probe[ip][iy+j][0] =  prr*trans[iyt+jt][0]
                                              - pri*trans[iyt+jt][1];  /* real */
                        probe[ip][iy+j][1] =  prr*trans[iyt+jt][1]
                                              + pri*trans[iyt+jt][0];  /* imag */
                    } /* end for(iy...) */
                }  /* end for(ix...) */
                fftwf_execute_dft( planPf, probe[ip], probe[ip] );
           }
    
            /*  multiplied by the propagator function */
            propagateW( probe[ip], propxr, propxi, propyr, propyi,
                kxp2, kyp2, (float)k2maxp, nxprobe, nyprobe );

        }  /* end  for( ip=... */

        /*  if this is a good thickness level then save the ADF or confocal signals
           - remember that the last level may be off by one layer with
              thermal displacements so special case it
        */
       
        /*  look at all values because they may not be in order */
        for( it = 0; it<nThick; it++ ) 
        if( fabs(ThickSave[it]-zslice)<fabs(0.5*deltaz)) {
            
            printf( "save ADF/confocal signals, thickness level %d\n", it);  /* diagnostic */
    
            /*  loop over all probes again */
#pragma omp parallel for private(ix,iy,j,idetect,prr,pri,delta,k2,cpix,phi,chi0,hr,hi,sum0,sum1,rx2,r2)
            for( ip=0; ip<npos; ip++) {
                /*  zero sum count */
                sum[ip] = 0.0;
                for(ix=0; ix<ndetect; ix++) detect[it][ix][ip] = 0.0;
    
                /*  sum intensity incident on the ADF detector
                        and calculate total integrated intensity
                    - changed detector limits to >= min and < max
                       so many concentric ADF detectors sum correctly
               7-jul-2011
                */
                for( ix=0; ix<nxprobe; ix++) {
                    j = ix*nyprobe;
                    for( iy=0; iy<nyprobe; iy++) {
                        prr = probe[ip][iy + j][0];
                        pri = probe[ip][iy + j][1];
                        delta = prr*prr + pri*pri;
                        sum[ip] += delta;
                        k2 = kxp2[ix] + kyp2[iy];
                        for( idetect=0; idetect<ndetect; idetect++) {
                            if( ADF == collectorMode[idetect] ) {
                                if( (k2 >= k2min[idetect] ) &&
                                (k2 < k2max[idetect] ) )
                                detect[it][idetect][ip] += delta;
                            }
                        }
                    } /* end for(iy..) */
                }  /* end for(ix...) */
                
                /*  transform back if confocal needed 
                    - use copy of probe so original can continue in use  */
                if( doConfocal == TRUE ) {
                    /*  allocate/deallocate here so openMP will work 
                        otherwise have to allocate nxprobe cpix arrays 
                        - a littel slow but a lot less memory */
                    cpix = (fftwf_complex*) fftwf_malloc( nxprobe*nyprobe * sizeof(fftwf_complex) );
                    if( NULL == cpix ) {
                        printf("Cannot allocate cpix array\n");
                        exit( EXIT_FAILURE );
                    }
                    sum0 = 0;
                    for( ix=0; ix<nxprobe; ix++) {
                        j = ix*nyprobe;
                        for( iy=0; iy<nyprobe; iy++) {
                            k2 = kxp2[ix] + kyp2[iy];
                            if( (k2 >= k2maxaC) && (k2 <= k2maxbC) ) {
                                phi = atan2( ky[iy], kx[ix] );
                                /*  offset defocus by zslice so both lens referenced to 
                                   entrance surface of specimen  */
                                chi0 = chi1*k2* ( (chi2C + chi3C*k2)*k2 - dfC + zslice
                                    + dfa2C*sin( 2.0*(phi-dfa2phiC) ) 
                                    + 2.0F*dfa3C*wavlen*sqrt(k2)*
                                    sin( 3.0*(phi-dfa3phiC) )/3.0 );
                                chi0= - chi0;
                                hr = (float) cos( chi0 );
                                hi = (float) sin( chi0 );
                                prr = probe[ip][iy + j][0];  /* real */
                                pri = probe[ip][iy + j][1];  /* imag */
                                cpix[iy + j][0] = prr*hr -pri*hi;
                                cpix[iy + j][1] = prr*hi +pri*hr;
                                sum0 += prr*prr + pri*pri;
                            } else {
                                cpix[iy + j][0] = 0.0F;
                                cpix[iy + j][1] = 0.0F;
                            }
                        }  /*  end for( iy... )  */
                    }  /*  end for( ix... )  */
            
                    fftwf_execute_dft( planPi, cpix, cpix );
                    scaleW( cpix, nxprobe, nyprobe );
                    
                    /* find normalization constant 
                      i.e. correct for constants in the FFT */
                    sum1 = 0.0;
                    for( ix=0; ix<nxprobe; ix++) {
                        j = ix*nyprobe;
                        for( iy=0; iy<nyprobe; iy++) {
                            prr = cpix[iy + j][0];
                            pri = cpix[iy + j][1];
                            sum1 += prr*prr + pri*pri;
                        }
                   }
            
                   /* integrate over real space detector and normalize */
                   for( ix=0; ix<nxprobe; ix++) {
                           rx2 = xoff[ip] - xp[ix];
                           rx2 = rx2*rx2;
                           j = ix*nyprobe;
                           for( iy=0; iy<nyprobe; iy++) {
                               r2 = yoff[ip] - yp[iy];
                               r2 = rx2 + r2*r2;
                               prr = cpix[iy + j][0];
                               pri = cpix[iy + j][1];
                               delta = prr*prr + pri*pri;
                               for( idetect=0; idetect<ndetect; idetect++) {
                                   if( CONFOCAL == collectorMode[idetect] ) {
                                     if( (r2 >= k2min[idetect] ) &&
                                               (r2 < k2max[idetect] ) )
                                         detect[it][idetect][ip] += delta*(sum0/sum1);
                                   }
                               }
                           }  /* end for( iy... )*/
                   }  /*  end for( ix....) */

                   fftwf_free( cpix );
                }  /* end if( doConfocal==TRUE) */
            
            }  /* end for( ip.. */

        }  /* end if( ((it...*/

       nslice++;
       zslice += deltaz;
       istart += na;

    }  /* end while( istart...) */

    free( ixoff );
    free( iyoff );
    free( xoff );
    free( yoff );

    return;

}/* end STEMsignals() */


/*--------------------- trlayer() -----------------------*/
/*  same subroutine in autoslic.c and autostem.c

  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  transr  = 2D array to get real part of specimen
        transmission function
  transi  = 2D array to get imag part of specimen
        transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

*/
void trlayer( const float x[], const float y[], const float occ[],
    const int Znum[], const int natom, 
    const float ax, const float by, const float kev,
    fftwf_complex *trans, const long nx, const long ny,
    double *phirms, long *nbeams, const float k2max )
{
    int idx, idy, i, j, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rmin2, sum, scale, scalex, scaley;

    const double rmax=3.0, rmax2=rmax*rmax; /* max atomic radius in Angstroms */

    scale = sigma( kev ) / 1000.0;  /* in 1/(volt-Angstroms) */

    scalex = ax/nx;
    scaley = by/ny;

    /* min radius to avoid  singularity */
    rmin = ax/((double)nx);
    r = by/((double)ny);
    rmin =  0.25 * sqrt( 0.5*(rmin*rmin + r*r) );
    rmin2 = rmin*rmin;

    idx = (int) ( nx*rmax/ax ) + 1;
    idy = (int) ( ny*rmax/by ) + 1;

    for( ix=0; ix<nx; ix++) {
        j = ix*ny;
        for( iy=0; iy<ny; iy++)
           trans[j++][0] = 0.0F;    /* real part trans[iy + ix*ny][0] */
    }
    
/*  paralleling this loop has little effect   */
/*#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,r,ixw,iyw,vz,rsq)*/
#pragma omp parallel for private(j,ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,ixw,iyw,vz,rsq)
    for( i=0; i<natom; i++) {
        ixo = (int) ( x[i]/scalex );
        iyo = (int) ( y[i]/scaley );
        nx1 = ixo - idx;
        nx2 = ixo + idx;
        ny1 = iyo - idy;
        ny2 = iyo + idy;

    /* add proj. atomic potential at a local region near its center
       taking advantage of small range of atomic potential */

        for( ix=nx1; ix<=nx2; ix++) {
            rx2 = x[i] - ((double)ix)*scalex;
            rx2 = rx2 * rx2;
            ixw = ix;
            while( ixw < 0 ) ixw = ixw + nx;
            ixw = ixw % nx;
            j = ixw*ny;
            for( iy=ny1; iy<=ny2; iy++) {
                rsq = y[i] - ((double)iy)*scaley;
                rsq = rx2 + rsq*rsq;
                if( rsq <= rmax2 ) {
                  iyw = iy;
                  while( iyw < 0 ) iyw = iyw + ny;
                  iyw = iyw % ny;
                  if( rsq < rmin2 ) rsq = rmin2;
                  /*r = sqrt( rx2 + r*r );
                  vz = occ[i] * vzatom( Znum[i], r ); slow */
                  vz = occ[i] * vzatomLUT( Znum[i], rsq );
                  trans[iyw + j][0] += (float) vz;
                }
            } /* end for(iy... */
       }  /* end for(ix... */

    } /* end for(i=0... */

    /* convert phase to a complex transmission function */
    sum = 0;
    for( ix=0; ix<nx; ix++) {
        j = ix*ny;
        for( iy=0; iy<ny; iy++) {
            vz = scale * trans[j][0];
            sum += vz;
            trans[j  ][0] = (float) cos( vz );  /* trans[iy+ix*ny][0] */
            trans[j++][1] = (float) sin( vz );
        }
    }

    *phirms = sum / ( ((double)nx)*((double)ny) );

    /* bandwidth limit the transmission function */
    *nbeams = 0;
    fftwf_execute_dft( planTf, trans, trans );
    for( ix=0; ix<nx; ix++) {
        j = ix*ny;
        for( iy=0; iy<ny; iy++) {
            k2 = ky2[iy] + kx2[ix];
            if (k2 < k2max) *nbeams += 1;
            else trans[iy + j][0] = trans[iy + j][1] = 0.0F;
        }
    }
    fftwf_execute_dft( planTi, trans, trans );
    scaleW( trans, nx, ny );
    
    return;

 }  /* end trlayer() */

/*------------------------- invert2D() ----------------------*/
/*
        rearrange pix with corners moved to center (a la FFT's)

         pix[ix][iy] = real array with image
         nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void invert2D( float** pix, long nx, long ny )
{
#define SWAP(a,b)       {t=a; a=b; b=t;}

        long ix, iy, ixmid, iymid;
        float t;

        ixmid = nx/2;
        iymid = ny/2;

        for( ix=0; ix<nx; ix++) 
        for( iy=0; iy<iymid; iy++)
                SWAP( pix[ix][iy], pix[ix][iy+iymid] );

        for( ix=0; ix<ixmid; ix++) 
        for( iy=0; iy<ny; iy++)
                SWAP( pix[ix][iy], pix[ix+ixmid][iy] );

#undef SWAP
}


