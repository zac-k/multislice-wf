/*
      *** atompot.c ***

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

  ANSI-C and TIFF version
  this version uses FFTW 3

  FFTW choses an optimum form of the FFT at run time so there
  is some variation in execution speed depending on what else 
  the CPU is doing during this planning stage

  see:   www.fftw.org

  on Windows file libfftw3f-3.dll must be in the PATH

  on Linux build as:
  g++ -O -fopenmp -o atompot atompot.cpp slicelib.o
                       floatTIFF.o -lfftw3f -lm

  Calculate the atomic potentials using tabulated electron
  scattering factors

  The units of the output are such that when multiplied by
  Lambda*(M/M0) the result is the phase shift of the transmitted
  electron wave (i.e. a simple sum of the unadulterated Born
  approx. scattering amplitudes).

  The calculation is performed in a 2D planar slice as a
  prelude to the multislice calculation of electron micrographs.

  Gaussian random displacements may be added to atom coordinate
  to simulate thermal vibrations (i.e. one snap shoot).

    rewritten in C  13-June-1995 ejk
    C version working 28-jun-1995 ejk
    fix bug in natom calculation 7-july-1995
    minor rearrangments 13-July-1995 ejk
    fix small bug in seval on end of range 16-july-1995 ejk
    converted to TIFF file format 8-apr-1996 ejk
    put fetab, featom() in slicelib.c 26-july-1996 ejk
    switch to new fe params 16-jan-1997 ejk
    switched to rangauss() random number generator 22-may-1997 ejk
    removed comma delimiters on input 5-july-1997 ejk
    remove option to set Fourier resolution and fix at max
        25-sep-1997 ejk
    switch to temperature input on thermal vibrations
        3-oct-1997 ejk
    cosmetic changed to compile on MAC/codeWarrior 5-mar-1998 ejk
    changed memory allocator functions 6-nov-1999 ejk
    change void main() to int main() for better portability
         21-jan-2000 ejk
    fix small input read error under new version 
     Linux (set iz=-1 before sscanf()) 22-jan-2000 ejk
    polished code in scamp() a little 8-feb-2006 ejk
    convert to GPL 3-jul-2008 ejk
    get return value of scanf() to remove warnings from gcc 4.4
      and convert to 4 char TAB size formatting 29-mar-2010 ejk
    convert to faster FFTW 4-apr-2010 ejk
    convert to floatTIFF.cpp 18-mar-2012 to 15-apr-2012 ejk

  This program is ANSI standard C++ and should be transportable
  although this is NOT guaranteed

  This source code is formatted for a tab size of 4.

*/



#include <cstdio>  /* ANSI C libraries */
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>

#include "../multislice/fftw3.h"         /* FFT routines from FFTW 3*/
//#include "floatTIFF.hpp"   /* file I/O routines in TIFF format */
#include "../multislice/floatTIFF.cpp"   /* file I/O routines in TIFF format */
//#include "slicelib.hpp"    /* misc. routines for multislice */
#include "../multislice/slicelib.cpp"    /* misc. routines for multislice */




const int NAMAX  = 2000;    /* max number of atoms of one type */
const int  NSMAX =  20;  /* max number of symmetry operations */
//const int  NZMAX =  103; /* max atomic number */
//const int  NCMAX =  132; /* characters per line to read */

int natom, nsx, nsy;  /* global vars for scamp() */
float *x, *y, *occ;

/* define subroutines defined at end of file */

void scamp( float, float, double*, double* );

int atomicPotential() 
{
    int ix,iy, nx, ny, ixmid, iymid, ns,
        i, j, nsym, jj, iz, is, lwobble;
    int ncellx, ncelly, ncellz;
    long ncoeff, ltime;
    unsigned long iseed;

    char filein[NCMAX], fileot[NCMAX], cline[NCMAX];

    float *pix, *kx, *ky, *wobble;
    float rmin, rmax, imin, imax, dx, dy, scalef;
    float *symx1, *symx2, *symy1, *symy2;
    double scampr, scampi, sum, runtime, fe, ky2, k2, k2max,
        scale, ax,by,cz,total1, total2, rx2, ry2,
        temperature, scalet;

    fftwf_complex *cpix;          /* complex fft */
    fftwf_plan  planTf;           /* FFTW plans */

    FILE *fp;

    floatTIFF myFile;

/*  the following are the chemical symbols for the periodic table */

    const char symbol[] = {
        " HHeLiBe B C N O FNeNaMgAlSi P SCl"
        "Ar KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBr"
        "KrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe"
        " IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTm"
        "YbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRn"
        "FrRaAcThPa UNpPuAmCmBkCfEsFmMdNoLr"
    };

    /*  Echo version date */
    printf("atompot version dated 15-apr-2012 EJK\n");
    printf("Copyright (C) 1998-2010 Earl J. Kirkland\n" );
    printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
        " under the GNU general public license\n\n" );

    printf( "calculate projected atomic potentials (to use in multislice)\n");
    printf( "    using FFTW\n\n");




	//************These operations have been removed because these variables will be set internally********

    /*  Get input file name etc. */
    /*printf("Name of file with input crystal data :\n");
    ns = scanf("%s", filein);

    printf("Name of file to get binary output"
          " of atomic potential :\n");
    ns = scanf("%s", fileot);

    printf("Real space dimensions in pixels Nx, Ny :\n");
    ns = scanf("%d %d", &nx, &ny);
	
    printf("Replicate unit cell by NCELLX,NCELLY,NCELLZ :\n");
    ns = scanf("%d %d %d", &ncellx, &ncelly, &ncellz);
    if( ncellx < 1 ) ncellx = 1;
    if( ncelly < 1 ) ncelly = 1;
    if( ncellz < 1 ) ncellz = 1;
*/

	/*Set pixel dimensions for multislice method*/
	nx = 256;
	ny = 256;
	
	/*Set replication factors*/
	ncellx = 10;
	ncelly = 10;
	ncellz = 10;





/*  read in parameters for Monte Carlo displacements
      to simulate thermal motion 
    get the initial seed from the time counter if it exists
    otherwise ask user for a seed 
   - this should get a different random number sequence each time the
      program is run 
*/
	/* **************************************************************************************************** 
       ************************Set lwobble = 1 if thermal displacement is to be used ********************** 
	   **************************************************************************************************** */
	//lwobble = askYN( "Do you want to add thermal displacements to atomic coord.?" );
    lwobble = 0;
	if( lwobble == 1 ) {
        printf("Temperature in degrees K:\n");
        ns = scanf("%lg", &temperature);
        ltime = (long) time( NULL );
        iseed = (unsigned) ltime;
        if( ltime == -1 ) {
            printf("Type initial seed for random number generator:\n");
            ns = scanf("%ld", &iseed);
        } else
            printf( "Random number seed initialized to %ld\n", iseed );
    }


/*  start timer */

    runtime = cputim();

/*  Start to read in specimen parameters
    remember to use fgets() in ReadLine() and NOT fscanf() 
    because fscanf() ignores newlines and 
        we need to sync on whole lines
*/

    /*
	fp = fopen( filein,"r" );
    if( fp == NULL ) {
        perror(" can't open crystal data input file");
        exit(0 );
    }
	
    ReadLine( fp, cline, NCMAX, filein );
    ns = sscanf( cline, "%lf %lf %lf", &ax, &by, &cz);
    printf("2D lattice constants= %f x %f Angstroms\n"
        " and propagation constant= %f Angstroms\n",
           ax, by, cz );*/

	/* Data for Si */

	ax = 3.8396; by = 6.4300; cz = 1.9198;//<----------------------------  Unit cell size

    if( (ncellx > 1) || (ncelly > 1) || (ncellz > 1 ) ) {
        ax = ax * ncellx;
        by = by * ncelly;
        cz = cz * ncellz;
        printf("Unit cell replicated to a= %g, b= %g, c=%g  Angstroms\n",
           ax, by, cz );
    }

    /* read in symmetry operations */

	/*
    ReadLine( fp, cline, NCMAX, filein );
    ns = sscanf( cline, "%d", &nsym);
	*/
	nsym = 0; //<-------------------------------------------------------  Number of symmetry operations
	/*
    if( nsym > 0 ) {
        symx1 = (float*) malloc1D( nsym, sizeof(float), "symx1" );
        symx2 = (float*) malloc1D( nsym, sizeof(float), "symx2" );
        symy1 = (float*) malloc1D( nsym, sizeof(float), "symy1" );
        symy2 = (float*) malloc1D( nsym, sizeof(float), "symy2" );
        for( i=0; i<nsym; i++)  {
            ReadLine( fp, cline, NCMAX, filein );
            ns = sscanf( cline, "%f %f %f %f",
                &symx1[i], &symx2[i], &symy1[i], &symy2[i]);
       }
    }
	/*

/*  Calculate misc constants (Nx*Ny added because FFT2D
    performs scaling) also adjust k2max for circular symmetry
    remember that k2max is still in units of Angstroms
*/
    k2max = 2.0*ax/nx;
    scale = 2.0*by/ny;
    if( scale > k2max ) k2max =scale;

    printf("Maximum symmetrical resolution set to %g Angstroms\n",
           k2max );
    k2max = 1.0 / (k2max * k2max);

    rx2 = (1.0/ax);  rx2 = rx2*rx2;
    ry2 = (1.0/by);  ry2 = ry2*ry2;
    scale = ( ((double)nx) * ((double)ny) ) /(ax*by);
    ixmid = nx/2;
    iymid = ny/2;

    ky = (float*) malloc1D( ny , sizeof(float), "ky" );
    kx = (float*) malloc1D( nx , sizeof(float), "kx" );

    /*  will only need half of these but doesn't take much CPU */
    for( iy=0; iy<ny; iy++) {
        if( iy < iymid )  ky[iy] = (float) iy;
        else ky[iy] = (float)(iy-ny);
    }

    for( ix=0; ix<nx; ix++) {
        if( ix < ixmid )  kx[ix] = (float) ix;
        else kx[ix] = (float)(ix-nx);
    }

    /*  allocate arays in FFTW format */
    pix = (float*) fftwf_malloc( nx*ny * sizeof(float) );
    cpix = (fftwf_complex*) fftwf_malloc( nx*(ny/2+1) * sizeof(fftwf_complex) );
    if( (NULL == pix) || (NULL == cpix) ) {
        printf("Cannot allocate memory in atompot\n" );
        exit( 0 );
    }
    /*  remember, FFTW uses the opposite sign convention */
    planTf = fftwf_plan_dft_c2r_2d( nx, ny, cpix, pix, FFTW_ESTIMATE );

    for( ix=0; ix<nx; ix++) 
    for( iy=0; iy<=iymid;iy++)  {
        cpix[ iy + ix*(iymid+1)][0] = 0.0F;  /* real */
        cpix[ iy + ix*(iymid+1)][1] = 0.0F;  /* imag */
    }

    /* ---- alloc misc arrays ----- */
    x      = (float*) malloc1D( NAMAX, sizeof(float), "x" );
    y      = (float*) malloc1D( NAMAX, sizeof(float), "y" );
    occ    = (float*) malloc1D( NAMAX, sizeof(float), "occ" );
    wobble = (float*) malloc1D( NAMAX, sizeof(float), "wobble" );

    /* ------ read in actual coordinates ------ */

    total2 = 0.0;
    printf("\n");
More:
        ReadLine( fp, cline, NCMAX, filein );
        iz = -1;    /* reset because some machines leave last value */
        ns = sscanf( cline, "%d", &iz);
        if( (strlen(cline) > 1) && (iz >= 1) && (iz <= NZMAX) ) {
            j = 0;
            total1 = 0.0;
            while( ReadLine( fp, cline, NCMAX, filein ) > 2 ) {
                ns = sscanf( cline, "%f %f %f %f", &occ[j], &x[j], &y[j], &wobble[j] );
                total1 = total1 +  occ[j] * (nsym+1) * ncellx * ncelly;
                if( nsym > 0) {
                    if( (j+nsym+1) > NAMAX) {
                        printf("Too many atoms\n");
                        printf("  Maximum allowed = %d\n",NAMAX);
                        fclose( fp );
                        exit(0);

                    }
                    for( jj=0; jj<nsym; jj++) {
                        occ[j+jj+1] = occ[j];
                        wobble[j+jj+1] = wobble[j];
                        x[j+jj+1] = (symx1[jj]*x[j] +symx2[jj] )/ncellx;
                        y[j+jj+1] = (symy1[jj]*y[j] +symy2[jj] )/ncelly;
                    }
                }
 
                x[j] = x[j] /ncellx;
                y[j] = y[j] /ncelly;
                j = j + nsym + 1;
                if( j > NAMAX ) {
                    printf("Too many atoms\n");
                    printf("  Maximum allowed = %d\n",NAMAX);
                    fclose( fp );
                    exit( 0 );
                }
            }  /* end while( ReadLine()>2 */

            natom = j ;
            is = 2*(iz-1);
            printf("%8.2f atoms with Z= %3d (%c%c)\n",
              total1, iz, symbol[is], symbol[is+1]);

            /*   if thermal vibrations are requested then we must expand
                ncellx,y here otherwise factor it inside scamp() 
                (its faster)
                 add random displacements to all positions
           */
            if( lwobble == 1 ) {
                if( (j-1)*ncellx*ncelly > NAMAX) {
                    printf("Too many atoms\n");
                    printf("  maximum allowed = %d\n",NAMAX);
                    exit( 0 );
                }

                if( ncellx > 1) {
                    for( ix=1; ix<ncellx; ix++)
                    for( i=0; i<natom; i++) {
                        x[j] = x[i] + ((float)ix)/((float)ncellx);
                        y[j] = y[i];
                        wobble[j] = wobble[i];
                        occ[j++] = occ[i];
                    }
                    natom = j - 1;
                }

                if( ncelly > 1) {
                    for( iy=1; iy<ncelly; iy++)
                    for( i=0; i<natom; i++) {
                        x[j] = x[i];
                        y[j] = y[i] + ((float)iy)/((float)ncelly);
                        wobble[j] = wobble[i];
                        occ[j++] = occ[i];
                    }
                    natom = j - 1;
                }

                /* scale thermal displacements to 300 deg. K 
                     and normalize to 3D - integrating over the z direction */
                scalet = sqrt( temperature / 300.0 );
                for( i=0; i<natom; i++) {
                    x[i] = x[i] + 
                        (float) (wobble[i] * scalet * rangauss( &iseed ) /ax);
                    y[i] = y[i] + 
                        (float) (wobble[i] * scalet * rangauss( &iseed ) /by);
                }

                nsx = 1;
                nsy = 1;

            } else {
                nsx = ncellx;
                nsy = ncelly;

            }  /* end if( lwobble == TRUE ) */

/*  calculate scattering amplitude in upper half plane (FFTW)
   NOTE zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for FFT, high freq is in the center

    - add extra complex conjugation to correct for the opposite
    sign convention used by FFTW

*/
            ncoeff = 0;
            for( iy=0; iy<=iymid; iy++) {
                ky2 = ky[iy]*ky[iy] * ry2;
                for( ix=0; ix<nx; ix++) {
                    k2 = kx[ix]*kx[ix] * rx2 + ky2;
                    if( k2 <= k2max) {
                        fe = scale * featom( iz, k2 );
                        scamp( kx[ix], ky[iy], &scampr, &scampi ) ;
                        cpix[ iy + ix*(iymid+1)][0] += (float) (scampr * fe);  /* real */
                        cpix[ iy + ix*(iymid+1)][1] += (float) (-scampi * fe); /* imag */
                        ncoeff++;
                    }

                }  /* end for(ix */
            } /* end for(iy */

            total2 = total2 + total1;
            goto More;
            
        }  /* end top if( iz>= 1) */
        
    /*  end loop over types */

    printf("\n   for a grand total of %10.2f atoms\n", total2);
    fclose( fp );

    fftwf_execute( planTf );  /*    inverse fft */

    /*  output results and find min and max to echo */

    if( nsym > 0 ) {
        free( symx1 );
        free( symx2 );
        free( symy1 );
        free( symy2 );
    }

    /*  copy to floatTIFF pix to output with */
    myFile.resize( nx, ny );
    myFile.setnpix( 1 );
    scalef = 1.0F/( nx*ny );  /*  scale factor missing in FFTW */
    rmin = rmax = scalef * pix[0];
    sum = 0.0;
    for( iy=0; iy<ny; iy++) {
        for( ix=0; ix<nx; ix++) {
            myFile(ix,iy)  =  scalef * pix[ iy + ix*ny];
            sum += myFile(ix,iy);
            if( myFile(ix,iy) < rmin ) rmin = myFile(ix,iy);
            if( myFile(ix,iy) > rmax ) rmax = myFile(ix,iy);
        }
    }

    myFile.zeroParam();
    myFile.setParam( pRMAX, rmax);
    myFile.setParam( pRMIN,  rmin);
    myFile.setParam( pC, (float) cz);
    myFile.setParam( pRES, (float) (1.0/sqrt(k2max)) );
    myFile.setParam( pDX,  dx = (float) ax / ((float)nx) );
    myFile.setParam( pDY,  dy = (float) by / ((float)ny) );

    imin = imax = 0;
    if( myFile.write( fileot, rmin,rmax, imin,imax, dx, dy ) != 1 )
        printf( "atompot cannot write an output file.\n");

    printf(" pix range %g to %g\n", rmin, rmax);
    printf("%ld fourier coeff. calculated in right half plane\n",
           ncoeff);
    printf("The average real space value was %g\n",
           sum /(((float)nx)*((float)ny)));
          
    printf("CPU time (excluding set-up) = %f sec.\n", cputim()-runtime );

    return EXIT_SUCCESS;

}  /* end main() */

/*--------------------- scamp() -----------------------*/
/*
  Complex specimen scattering amplitude

  kx,ky = real scattering vectors
  scampr, scampi = real and imag. parts of scattering amplitude
           returned by this routine

  global variables used:
  x[natom], y[natom] = real array with atomic coord.
  occ[natom] = real array with occupations at each x,y
  natom      = integer number of atoms (all with same Z)
  nsx, nsy = integer number of cells to replicate in x,y

*/
void scamp( float kx, float ky, double *scampr, double *scampi )
{
    int i, ixc, iyc;
    double pi2, w1, w;
    double  scalex, scaley;

    double phasexr, phasexi, phaseyr, phaseyi;

    pi2 = 8.0 * atan( 1.0 );
    
    scalex = pi2*kx / ((double)nsx);
    scaley = pi2*ky / ((double)nsy);

    *scampr = *scampi = 0.0;

/*  use a trick to sum nsx + nsy terms
    instead of nsx * nsy
    - cos()/sin() could be done recursive for speed
*/
    for( i=0; i<natom; i++) {

        phasexr = phasexi = 0.0;
        w1 = pi2 * kx * x[i];
        for( ixc=0; ixc<nsx; ixc++) {
            w = w1 + ixc*scalex;
            phasexr += cos(w);
            phasexi += sin(w);
        }

        phaseyr = phaseyi = 0.0;
        w1 = pi2 * ky * y[i];
        for( iyc=0; iyc<nsy; iyc++) {
            w = w1 + iyc*scaley;
            phaseyr += cos(w);
            phaseyi += sin(w);
        }

        *scampr += ( phasexr*phaseyr - phasexi*phaseyi ) * occ[i];
        *scampi += ( phasexr*phaseyi + phasexi*phaseyr ) * occ[i];

    } /* end for(i=0... */

    return;

 }  /* end scamp() */
