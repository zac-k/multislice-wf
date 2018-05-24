/*      *** probe.cpp ***

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

------------------------------------------------------------------------

    ANSI-C version
    this version uses FFTW 3

    see:   www.fftw.org

    on Windows file libfftw3f-3.dll must be in the PATH

    on Linux build as:
    g++ -O -o probe probe.cpp slicelib.o
                       floatTIFF.cpp -lfftw3f

    Calculate a focused probe wavefunction in real space

    this file is formatted for a tab size of 4 characters

    rewritten in C 6-dec-1995 ejk
    fixed sign error in aberration function 1-mar-1997 ejk
    removed commas from keyboard input 3-oct-1997 ejk
    updated memory allocation routines 20-nov-1999 ejk
    change void main() to int main() for better portability
         22-jan-2000 ejk
    add Cs5=5th order aberration and astigmatism  19-jul-2005 ejk
    small cosmetic changes 18-jul-2007 ejk
    convert to GPL 3-jul-2008 ejk
        convert to large list aber. with coma option 23-nov-2008 ejk
    get return value of scanf() to remove warnings from gcc 4.4
      and convert to 4 char TAB size formatting 10-apr-2010 ejk
    convert to FFTW 9-may-2010 ejk
    fix C34a,b terms 10-may-2010 ejk
    change astig. parameters to a,b from mag.+angle 30-may-2010 ejk
    add more aberations to 5th order 30-jun-2010 to 4-jul-2010 ejk
    add probe size calculation 5-jul-2010 ejk
    split up into subroutine 16-jul-2010 ejk
    fix a few things in prbSize() 5-sep-2010 ejk
    switch to storing multipole aberr. in param[] to save in
       output file 2-may-2011 ejk
    add multiMode in chi() to avoid extra calculations if there
       are no multipole aberrations 8-may-2011 ejk
    convert to C++ and floatTIFF.cpp  22-mar-2012 ejk

*/

#include <cstdio>  /*  ANSI-C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include "fftw3.h"         /* FFT routines from FFTW 3*/
#include "slicelib.hpp"    /* define parameter offsets */
#include "floatTIFF.hpp"   /* file I/O libraries */

#define MANY_ABERR      /*  define to include many aberrations */

const int NCMAX=   1024; /* characters per line */

/* global data */
float *kx, *ky, *xpos, *ypos, *kx2, *ky2;
fftwf_plan planTi;    /* FFTW plans */

/*  define subroutines at bottom of file */
int makeProbe( fftwf_complex *cpix, int nx, int ny,
          double xp, double yp, float p[], 
          double wavlen, double k2max, double pixel, int multiMode,
          int ismoth );
double prbSize( float** pixsq, int nx, int ny,
    double xp, double yp, double ax, double by );


/*------------------------ main() ---------------------*/
int main()
{
    char fileout[NCMAX], cline[NCMAX];
    const char version[] = "22-mar-2012 (ejk)";
    int ix, iy, nx, ny, ixmid, iymid, i, ismoth, npixels, ns,
        done, status, multiMode, NPARAM;
    float rmin, rmax, aimin, aimax, p2;
    float *param, pixr, pixi, **pixsq;
    double k2max, keV, wavlen, ax, by, rx, ry,
        rx2, ry2, pi, dx, dy, scale, scale2, pixel,
        Cs3, Cs5, df, sum, time, apert;
    double  x;

    fftwf_complex *cpix;  /* complex arrays for fft */

    floatTIFF myFile;

/*  Echo version date etc.  */
    printf( "probe version dated %s\n", version );
    printf("Copyright (C) 1998-2011 Earl J. Kirkland\n" );
    printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
        " under the GNU general public license\n\n" );

#ifdef MANY_ABERR
    printf( "calculate a focused probe wave function including multiple aberr.\n\n");
#else
    printf( "calculate a focused probe wave function\n\n");
#endif

    pi = 4.0 * atan( 1.0 );

    /*  memory to store parameters */
    NPARAM = myFile.maxParam();
    param = (float*) malloc1D( NPARAM, sizeof(float), "probe-param" );
    for( i=0; i<NPARAM; i++) param[i] = 0.0F;

/* ---- Get desired image size, parameters etc. ------------- */

    printf("Name of file to get focused probe wave function:\n");
    ns = scanf("%s", fileout );

    printf("Desired size of output image in pixels Nx,Ny:\n");
    ns = scanf("%d %d", &nx, &ny );

    printf("Size of output image in Angstroms ax,by:\n");
    ns = scanf("%lf %lf", &ax, &by );

    printf("Probe parameters, V0(kv), Cs3(mm), Cs5(mm),"
           " df(Angstroms), apert(mrad):\n");
    ns = scanf("%lg %lg %lg %lg %lg",
          &keV, &Cs3, &Cs5, &df, &apert );
    param[pDEFOCUS] = (float) df;
    param[pCS]  = (float) ( Cs3*1.0e7 );
    param[pCS5] = (float) ( Cs5*1.0e7 );

    printf("Type 1 for smooth aperture:\n");
    ns = scanf("%d", &ismoth );

    printf("Probe position x,y in Ang.:\n");
    ns = scanf("%lf %lf", &dx, &dy );

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
            ns = scanf( "%lg", &x );
            /* printf("%s, %f\n", cline, x );  testing */
            status = readCnm( cline, param, x );        
            if( status < 0 ) {
                printf( "unrecognized aberration, exit...\n");
                exit( EXIT_SUCCESS );
            } else multiMode = 1;
        }
    } while( !done );

#endif

/* ------- Calculate misc constants ------------ */

    time = cputim( );
    
    rx  = 1.0/ax;
    rx2 = rx * rx;
    ry  = 1.0/by;
    ry2 = ry * ry;
    
    ixmid = nx/2;
    iymid = ny/2;
    
    wavlen = wavelength( keV );
    printf("electron wavelength = %g Angstroms\n", wavlen);

    k2max = apert*0.001/wavlen;
    k2max = k2max * k2max;

/* ------- allocate memory ------------ */

    pixsq = (float**) malloc2D( nx, ny, sizeof(float), "pixsq" );

    kx   = (float*) malloc1D( nx, sizeof(float), "kx" );
    kx2  = (float*) malloc1D( nx, sizeof(float), "kx2" );
    xpos = (float*) malloc1D( nx, sizeof(float), "xpos" );
    freqn( kx, kx2, xpos, nx, ax );

    ky   = (float*) malloc1D( ny, sizeof(float), "ky" );
    ky2  = (float*) malloc1D( ny, sizeof(float), "ky2" );
    ypos = (float*) malloc1D( ny, sizeof(float), "ypos" );
    freqn( ky, ky2, ypos, ny, by );


    /*------  make FFTW arrays and plans ------- */
    cpix = (fftwf_complex*) fftwf_malloc( nx*ny * sizeof(fftwf_complex) );
    if( NULL == cpix ) {
        printf("Cannot allocate wave array cpix\n");
        exit( EXIT_FAILURE );
    }
    /* remember FFTW has inverse sign convention */
    planTi = fftwf_plan_dft_2d( nx, ny, cpix, cpix, 
        FFTW_FORWARD, FFTW_ESTIMATE );  /* inverse in place */

    /* --------- calculate probe wavefunction -------- */
    pixel = ( rx2 + ry2 );
    npixels = makeProbe( cpix, nx, ny, dx, dy, 
        param, wavlen, k2max, pixel, multiMode, ismoth );
    printf("there were %d pixels inside the aperture\n", npixels );

    /* -----  copy back for output ----- */
    sum = 0.0;
    myFile.resize( 2*nx, ny);
    myFile.setnpix( 2 );
     for( ix=0; ix<nx; ix++)
    for( iy=0; iy<ny; iy++) {
        myFile(ix,iy)    = pixr = cpix[iy + ix*ny][0];  // real
        myFile(ix+nx,iy) = pixi = cpix[iy + ix*ny][1];  // imag
        pixsq[ix][iy] = p2 = pixr*pixr + pixi*pixi;
        sum +=  p2;
    }

/* ----- Normalize probe intensity to unity ------------ */

    scale = 1.0 / sum;
    scale = scale * ((double)nx) * ((double)ny);
    scale = (double) sqrt( scale );
    scale2 = scale*scale;

    for( ix=0; ix<nx; ix++) 
       for( iy=0; iy<ny; iy++) {
        myFile(ix,iy)    *= (float) scale;
        myFile(ix+nx,iy) *= (float) scale;
        pixsq[ix][iy]    *= (float) scale2;
    }

/*------- Output results and find min and max to echo ---------------
*/
    rmin  = myFile.min(0);    // real part
    rmax  = myFile.max(0);
    aimin = myFile.min(1);   // imaginary
    aimax = myFile.max(1);

    param[pRMAX] = rmax;
    param[pIMAX] = aimax;
    param[pRMIN] = rmin;
    param[pIMIN] = aimin;
    param[pDEFOCUS]= (float) df;
    param[pDX]= (float) (ax / nx);
    param[pDY]= (float) (by / ny);
    param[pENERGY]= (float) keV;
    param[pWAVEL]= (float) ( sqrt(k2max) * wavlen);
    param[pCS]= (float) Cs3;
    param[27]= (float) dx;
    param[28]= (float) dy;

    for( i=0; i<NPARAM; i++) myFile.setParam( i, param[i] );  // not very efficient

    if( myFile.write( fileout, rmin, rmax, aimin, aimax, (float) dx, (float) dy ) != 1 )
        printf( "probe cannot write an output file.\n");

    printf( "Pix range %15.7g to %15.7g real,\n"
        "      and %15.7g to %15.7g imaginary\n",
        rmin, rmax, aimin, aimax );

/*------- calculate probe size ---------------*/

    x = prbSize( pixsq, nx, ny, dx, dy, ax, by );
    printf("probe size (FWHM-II) = %g Ang.\n", x);

/*------- exit ---------------*/

    time = cputim() - time;
    printf("\nCPU time = %f sec\n", time );

    return EXIT_SUCCESS;

}  /* end main() */

/*   ------------------ makeProbe() ----------------------
    calculate aberration limited probe wavefunction

    NOTE zero freq. is in the bottom left corner and
    expands into all other corners - not in the center
    this is required for FFT

  input:
    cpix[] = fftw array to get wave function
    nx, ny = size of wavefunciton in pixels
    xp, yp = probe position in Ang.
    p[]    = aberration values
    k2max  = max k^2
    pixel  = smoothing range
    multiMode = if not 0 then include multipole aberrations
    ismoth = flag to control smoothing at the edge

  returned value = number of pixels in the aperture

  assumed globals:
      planeTi
      kx[], kx2[], xpos[], ky[], ky2[], ypos[]

*/

int makeProbe( fftwf_complex *cpix, int nx, int ny,
          double xp, double yp, float p[], 
          double wavlen, double k2max, double pixel, int multiMode,
          int ismoth )
{ 
    int ix, iy, npixels;
    double alx, aly, k2, chi0, pi, dx2p, dy2p;

    /*    PIXEL = diagonal width of pixel squared
        if a pixel is on the aperture boundary give it a weight
        of 1/2 otherwise 1 or 0
    */
    npixels = 0;
    pi = 4.0 * atan( 1.0 );

    dx2p = 2.0*pi*xp;
    dy2p = 2.0*pi*yp;

    for( iy=0; iy<ny; iy++) {
        aly = wavlen * ky[iy];  /* y component of angle alpha */
        for( ix=0; ix<nx; ix++) {
            k2 = kx2[ix] + ky2[iy];
            alx = wavlen * kx[ix];  /* x component of angle alpha */
            if ( ( ismoth != 0) && 
                ( fabs(k2-k2max) <= pixel) ) {
                   chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode )
                           - ( (dx2p*kx[ix]) + (dy2p*ky[iy]) );
                cpix[iy + ix*ny][0] = (float) ( 0.5 * cos(chi0));  /* real */
                cpix[iy + ix*ny][1] = (float) (-0.5 * sin(chi0));  /* imag */
                /* printf("smooth by 0.5 at ix=%d, iy=%d\n", ix, iy ); */
            } else if ( k2 <= k2max ) {
                   chi0 = (2.0*pi/wavlen) * chi( p, alx, aly, multiMode )
                           - ( (dx2p*kx[ix]) + (dy2p*ky[iy]) );
                cpix[iy + ix*ny][0] = (float)  cos(chi0);  /* real */
                cpix[iy + ix*ny][1] = (float) -sin(chi0);  /* imag */
                npixels++;
            } else {
                cpix[iy + ix*ny][0] = cpix[iy + ix*ny][1] = 0.0F;
            }
        }  /* end for( ix=0... )  */
    }  /* end for( iy=0... )  */


/* ------- inverse transform back to real space ------------ */

    fftwf_execute_dft( planTi, cpix, cpix );
    scaleW( cpix, nx, ny );

    return( npixels );

}  /* end makeProbe()  */

/* -------------------  prbsSize() -------------------
   calculate probe size 

  input:
     pixsq[][] = intensity in probe
     nx, ny = size of probe in pixels
     xp, yp = original probe position (in Ang.)
     ax, by = size of cell in Ang.

  return value:  size of probe in Ang.

*/
double prbSize( float** pixsq, int nx, int ny,
    double xp, double yp, double ax, double by )
{
    int ix, iy, ir, nr;
    double *ivsr, dr, scale, rx, ry, x, y, y2;

    dr = ax/nx;
    scale = by/ny;
    if( scale < dr ) dr = scale;  /*  smallest radial spacing */
    dr = 0.5*dr;                  /*  use sub pixel sampling */
    if( ny > nx ) nr = ny;
    else nr = nx;
    nr = 2*nr;
    ivsr = (double*) malloc1D( nr, sizeof(double), "ivsr");
    for( ir=0; ir<nr; ir++) ivsr[ir] = 0.0;

    /*  find curent vs. radius = azimuthal average */
    /*    orig. probe position = (dx,dy)  */
    rx = ax/nx;
    ry = by/ny;
    for( iy=0; iy<ny; iy++) {
        y = (iy*ry) - yp;
        y2 = y*y;
        for( ix=0; ix<nx; ix++) {
            scale = pixsq[ix][iy];
            x = (ix*rx) - xp;
            ir = (int)( sqrt( x*x + y2 )/dr + 0.5);
        if( ir < 0 ) ir = 0;
        if( ir >= nr ) ir = nr-1;
            ivsr[ir] += scale;
       }
    }

    /*  integrate current */
    for( ir=1; ir<nr; ir++) ivsr[ir] += ivsr[ir-1];

    /*  find half current radius */
    scale = 0.5 * ivsr[nr-1];
    for( ir=0; ir<nr; ir++) {
        ix = ir;
        if( ivsr[ir] > scale ) break;
    }

    y = ivsr[ix] - ivsr[ix-1];   /* interpolate */
    if( fabs(y) > 0.0 ) y = (scale-ivsr[ix-1])/y;
    else y = 0.0;
    x = 2.0 * ( (ix-1)*dr + y*dr );

    free( ivsr );

    return( x );

}  /*  end prbSize() */
